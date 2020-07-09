class Plot:
    
    def __init__(self):
        pass

    def plot(self,filename):
        try:
            import matplotlib.pyplot as plt
        except:
            print("ImportError: matplotlib.pyplot @ plot")

        try:
            import numpy as np
        except:
            print("ImportError: Numpy @ plot")

        plt.rc('text', usetex=True)
        plt.style.use('seaborn-white')

        #-------------------------------------------------------------------------------------------------

        linestyle = {"densly_dashdotted":  (0, (3, 1, 1, 1)), }
        data = np.genfromtxt(filename)

        _, ax = plt.subplots(1, 2)
        ax[0].tick_params(labelsize=10, direction='in', length=5,
                    width=.5, colors='k', grid_color='r', grid_alpha=0.5)

        ax[0].plot(data[:, 0]/1000, data[:, 1], color="red", linestyle="-")
        ax[0].set_xlabel("$r~[10^3 km]$")
        ax[0].set_ylabel("$M/M_{\\odot}$")
        ax[0].set_title("$Mass-Distribution$")

        ax[1].tick_params(labelsize=10, direction='in', length=5,
                        width=.5, colors='k', grid_color='r', grid_alpha=0.5)
                        
        ax[1].plot(data[:, 0]/1000, data[:, 2], color="red",linestyle=linestyle["densly_dashdotted"])
        ax[1].set_xlabel("$r~[10^3 km]$")
        ax[1].set_ylabel("$P/\\epsilon_0$")
        ax[1].set_title("$Pressure-Distribution$")

        plt.show()
