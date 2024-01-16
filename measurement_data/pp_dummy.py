import numpy as np
import matplotlib.pyplot as plt
import sys
#sys.path.append('./RT')




class m_data:

    def __init__(self):
        self.time = []
        self.voltage = []
        self.current = []
        self.charge_cap = []
        self.discharge_cap = []
        self.temperature = []
        self.delta = 1e0          # [-delta, delta] for current in which battery is supposed to by idle


    def from_csv(self, fi):

        skipper = self.skip_first_lines(fi)

        values = np.genfromtxt(
            skipper, 
            dtype="S32, f16, f16, f16, f16, f16, f16,",
            missing_values="",
            delimiter=",",
            invalid_raise=False, 
            filling_values=0,
        )

        self.ID = str(values[0][0])


        for elem in values:
            self.time.append(elem[1])
            self.voltage.append(elem[2])
            self.current.append(elem[3])
            self.charge_cap.append(elem[4])
            self.discharge_cap.append(elem[5])
            self.temperature.append(elem[6])


    def from_myself(self, add_ID, i_start, i_end):
        new_self = m_data()

        new_self.ID = self.ID + str(add_ID)

        new_self.time = self.time[i_start:i_end]
        new_self.voltage = self.voltage[i_start:i_end]
        new_self.current = self.current[i_start:i_end]
        new_self.charge_cap = self.charge_cap[i_start:i_end]
        new_self.discharge_cap = self.discharge_cap[i_start:i_end]
        new_self.temperature = self.temperature[i_start:i_end]

        return new_self
    

    def skip_first_lines(self, fi):
        i = 1
        for line in fi:
            if i < 9:
                i += 1
                continue 
            yield line


    def my_subplot(self, ax, value, title):
        ax.set_title(title)
        ax.plot(self.time, value)
        ax.locator_params(nbins=3)
        ax.set_xlabel('time / s')


    def plot_me(self):
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=5, ncols=1)
        self.my_subplot(ax1, self.voltage, 'voltage / mV')
        self.my_subplot(ax2, self.current, 'current / mA')
        self.my_subplot(ax3, self.charge_cap, 'charge_cap / mAh')
        self.my_subplot(ax4, self.discharge_cap, 'discharge_cap / mAh')
        self.my_subplot(ax5, self.temperature, 'tamperature / ˚C')
        fig.suptitle("ID: " + self.ID)
        plt.show()

    def find_cycle_codes(self):
        i_start, state = [], []

        #state:
        #   0 = discharge
        #   1 = idle
        #   2 = charge
        #   3 = unknown

        state_tmp = 3
        start_tmp = 0
        d = self.delta
        time_out = 0

        for i in range(len(self.current)):
            if time_out < 1:
                if state_tmp == 0:
                    # find idle: 
                    if self.current[i] > -d:
                        state_tmp = 1
                        start_tmp = i

                elif state_tmp == 1:
                    # find discharge:
                    if self.current[i] < -d:
                        state_tmp = 0
                        state.append(0)
                        i_start.append(int((i + start_tmp) / 2))

                    # or find charge:
                    elif self.current[i] > d:
                        state_tmp = 2
                        state.append(2)
                        i_start.append(int((i + start_tmp) / 2))

                    else:
                        continue

                    time_out = 180

                elif state_tmp == 2:
                    # find idle: 
                    if self.current[i] < d:
                        state_tmp = 1
                        start_tmp = i

                elif state_tmp == 3:
                    # find idle: 
                    if self.current[i] < d:# and self.current[i] > -d:
                        state_tmp = 1

                else:
                    print("Error: variable state in m_data.find_cycle_codes in unknown state")
            else: 
                time_out -= 1

        return [i_start, state]
    

    def find_cycles(self):
        [i_start, state] = self.find_cycle_codes()
        charge_cycles = []
        discharge_cycles = []

        i_stop = i_start[1:]
        i_stop.append(len(self.current))
        i_stop = [elem - 1 for elem in i_stop]


        for i in range(len(i_start)):
            if state[i] == 0:
                discharge_cycles.append(self.from_myself(i, i_start[i], i_stop[i]))


            elif state[i] == 2:
                charge_cycles.append(self.from_myself(i, i_start[i], i_stop[i]))



        return [charge_cycles, discharge_cycles]
    
    def capacity(self):
        Power = np.zeros(len(self.current))
        for i in range(len(self.current)):
            Power[i] = abs(self.current[i]) * self.voltage[i]

        self.cap = np.trapz(Power, x = self.time) * 1e-9/3.6    #Wh
        return self.cap




class SOC_data:

    def __init__(self, cycles, total_cap):
        self.resol = 100
        self.ID = cycles.ID
        [self.SOC, self.Cap, self.time] = self.capacity_array(cycles, total_cap)
        [self.voltage, self.current, self.temperature] = self.values_at_Cap_time(cycles)


    def my_subplot(self, ax, value, title):
        ax.set_title(title)
        ax.plot(self.SOC, value)
        ax.locator_params(nbins=3)
        ax.set_xlabel('SOC')
        max_val = abs(max(value, key=abs))
        ax.axis([1.05, -0.05, min(value) - 0.1*max_val, max(value) + 0.1*max_val])



    def plot_me(self):
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
        plt.gca().invert_xaxis()
        self.my_subplot(ax1, self.voltage, 'voltage / mV')
        self.my_subplot(ax2, self.current, 'current / mA')
        self.my_subplot(ax3, self.Cap, 'Capacity / mA')
        self.my_subplot(ax4, self.temperature, 'tamperature / ˚C')
        fig.suptitle("ID: " + self.ID)

        plt.show()
        

    def values_at_Cap_time(self, cycle):

        voltage = np.zeros(shape=[self.resol, 1])
        current = np.zeros_like(voltage)
        temperature = np.zeros_like(voltage)
        idx = 0

        for i in range(self.resol):

            vol_tmp, cur_tmp, temp_tmp, count = 0, 0, 0, 0
            while(cycle.time[idx] < self.time[i] and idx < len(cycle.time)):
                vol_tmp += cycle.voltage[idx]
                cur_tmp += cycle.current[idx]
                temp_tmp += cycle.temperature[idx]
                count += 1
                idx += 1

            voltage[i] = vol_tmp / count
            current[i] = cur_tmp / count
            temperature[i] = temp_tmp / count

        return [voltage, current, temperature]
    

    def capacity_array(self, cycle, total_cap):

        
        Cap = np.zeros(shape=(self.resol, 1))
        Cap_time = np.zeros_like(Cap)
        SOC = np.zeros_like(Cap)

        Cap_tmp, idx = 0, 1
        P_old = abs(cycle.current[0]) * cycle.voltage[0]
        comparing_cap = total_cap*0.9999*1e9*3.6

        for i in range(self.resol):

            while(Cap_tmp < ((i+1)/self.resol*comparing_cap) and idx < len(cycle.current)):
                P = abs(cycle.current[idx]) * cycle.voltage[idx] 
                dt = cycle.time[idx] - cycle.time[idx-1]
                Cap_tmp += (P + P_old) / 2 * dt
                idx += 1
                P_old = P 

            Cap_time[i] = cycle.time[idx]
            SOC[i] = 1 - (i+1)/self.resol
            Cap[i] = Cap_tmp * 1e-9/3.6    #Wh

        return [SOC, Cap, Cap_time]


def capacity_investigation(charge_cycles, discharge_cycles):
    charge_cap = []
    for elem in charge_cycles:
        #elem.plot_me()
        charge_cap.append(elem.capacity())


    discharge_cap = []
    for elem in discharge_cycles:
        #elem.plot_me()
        discharge_cap.append(elem.capacity())

    return [charge_cap, discharge_cap]


def plot_field(x, y, z_field, name):
    plt.plot(y[:2], z_field[0, :2], 'o', label = (str(x[0]) + "°C"))
    plt.plot(y[:3], z_field[1, :3], 'o', label = (str(x[1]) + "°C"))
    plt.plot(y[:], z_field[2, :], 'o', label = (str(x[2]) + "°C"))
    plt.plot(y[:], z_field[3, :], 'o', label = (str(x[3]) + "°C"))
    plt.xlabel("C_rate / C")
    plt.ylabel("Capacity / Wh")
    plt.title(name + " over C_rate")
    plt.legend()
    
    #plt.axis([0, 11,0, 0.5])       # for relative loss
    plt.show()


def print_data(x, y, data, name):
    print("x-values:")
    print(x)
    print("y_values:")
    print(y)
    print(name + " - data:")
    print(data)




### main ###

## find battery cycles   
csv = open("dummy_battery_data.csv", "rb")
battery = m_data()
battery.from_csv(csv)
battery.plot_me()

[battery_charge_cycles, battery_discharge_cycles] = battery.find_cycles() 
for elem in battery_charge_cycles:
    elem.plot_me()
for elem in battery_discharge_cycles:
    elem.plot_me()

[battery_charge_cap, battery_discharge_cap] = capacity_investigation(battery_charge_cycles, battery_discharge_cycles)


#battery_SOC_data = SOC_data(battery_charge_cycles, battery_discharge_cycles)

battery_SOC_data = []
for elem in range(1):
    tmp = SOC_data(battery_charge_cycles[elem], battery_discharge_cycles[elem])
    battery_SOC_data.append(tmp)

### end main ###
