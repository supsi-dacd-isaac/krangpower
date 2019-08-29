import krangpower as kp
import numpy as np


# In this example, we will simulate a battery with a toy logic that decides what to do according to its state and
# a price signal.


# -------------------------------------
# CONSTANTS
# -------------------------------------
PRICEPERIOD = 48
um = kp.UM


# -------------------------------------
# PRICE DEFINITION
# -------------------------------------
# we define a simple sinusoidal price with 12 hour period.
def price(t_hour):
    return np.sin(t_hour/PRICEPERIOD/2/np.pi) + np.random.uniform(-0.7, 0.7)


# -------------------------------------
# BATTERY SMART LOGIC DEFINITION
# -------------------------------------
# we have to define the logic we want for the battery into a DecisionModel inherited class.
class BatteryManager(kp.DecisionModel):

    def __init__(self, capacity=10*um.kWh, init_soc=5*um.kWh, pwr=3*um.kW, aggressiveness=2):
        self.avg_price = 0
        self.aggro = aggressiveness
        self.nsamples = 0
        self.capacity = capacity
        self.SOC = init_soc
        self.pwr_limit = pwr
        self.money_gained = 0
        self.pf = 0.985

    def decide_pq(self, oek, mynode):  # you have to override "decide_pq".
        # the logic implemented here decides power according to how low is the price with respect to the recorded
        # average and how much the battery is far from being filled at half capacity

        # get the price
        hour = oek.brain.Solution.DblHour()  # returns the simulated hour
        el_price = price(hour.magnitude)  # this is the sell/buy price

        # update the average
        self.avg_price = (self.avg_price * self.nsamples + el_price) / (self.nsamples + 1)
        self.nsamples += 1

        # calculate conveniency and hunger (both are in [-1,1])
        conveniency = self.avg_price - el_price
        hunger = (self.capacity/2 - self.SOC)/(self.capacity/2)

        # calculate the power as the sum of conveniency and hunger, clipped to the battery power limits. The logic
        # operates in GENERATOR convention: positive power is produced, negative power is absorbed.
        newpower = np.clip(
            - (conveniency + hunger) * self.aggro,
            - self.pwr_limit.magnitude,
            self.pwr_limit.magnitude
        )

        energy_exchanged = newpower * um.kW * oek.get('stepsize')['stepsize']

        # don't forget to update your SOC!
        self.SOC += -energy_exchanged  # negative energy exchanged CHARGES the battery
        self.SOC = np.clip(self.SOC, 0*um.kWh, self.capacity)

        # we update the gain/expense
        self.money_gained += (energy_exchanged.to('kWh').magnitude * el_price)

        # you must give back P, Q
        return newpower * um.kW, 0.01 * newpower * um.kW


# -------------------------------------
# CREATING THE CIRCUIT
# -------------------------------------
# creating an appropriate voltage source. We want to simulate a 400V, 50Hz grid.
vs = kp.Vsource(basekv=0.4, frequency=50.0, Isc1=5)

# instantiating the circuit
circuit = kp.Krang('market_example', voltage_source=vs)

# we add a couple lines
circuit['sourcebus', 'bus_1'] << kp.Line(length=50 * um.m).aka('line_1')
circuit['bus_1', 'bus_2'] << kp.Line(length=50 * um.m).aka('line_2')

# definition and addition of our "smart battery" to bus_2
mylogic = BatteryManager()
circuit['bus_2', ] << kp.FourQ(kV=0.4).aka('controlled_battery') * mylogic

# setting the simulation step at 1 hr, tot sim time at three days
circuit.set(stepsize=1 * um.hr, number=72)


# -------------------------------------
# DEFINING WHAT RESULTS TO RECORD
# -------------------------------------
# in order to get interesting results, you have to define functions that take the Krang as argument and calculate
# interesting stuff. You can also refer to other objects in the namespace.

def my_power(oek):
    pwr = sum(oek['controlled_battery'].Powers()[0][0:3]).magnitude  
    # we have to call magnitude, because they're all
    # Pint Quantities!
    return np.real(pwr)


def my_soc(_oek):
    return (mylogic.SOC / mylogic.capacity).magnitude


def voltage_at_main(oek):
    return (np.sum(np.abs((oek['sourcebus'].Voltages()[0:3]))) / 3).magnitude


# -------------------------------------
# SOLVE AND DISPLAY
# -------------------------------------
# evalsolve solves the prescribed steps of simulation and at each step evaluates the functions you pass to it.
pwr_hist = circuit.evalsolve(my_power, my_soc, voltage_at_main, as_df=False)

# we get a dict whose keys are the function names and whose values are lists of the results returned by the functions.
# print(pwr_hist)

# we display how many money units we gained with our wise battery operation.
print('\nBalance:')
print(mylogic.money_gained)
