class Plot:
    
    def __init__(self):
        pass
    
    #-------------------------------------------------------------------------------------------------

    def double_plot(self,filename):
        try:
            import matplotlib.pyplot as plt
        except:
            print("ImportError: matplotlib.pyplot @ plot")

        try:
            import numpy as np
        except:
            print("ImportError: Numpy @ plot")

        #-------------------------------------------------------------------------------------------------

        plt.rc('text', usetex=True)
        plt.style.use('seaborn-white')

        #-------------------------------------------------------------------------------------------------

        lw = 2
        color = {
            "indigo": "#4B0082",
            "royal_blue": "#4169E1",
            "crimson": "#B22222",
            "salmon": "#FA8072",
        }

        #-------------------------------------------------------------------------------------------------

        data = np.genfromtxt(filename)

        _, ax = plt.subplots(1, 2)
        ax[0].tick_params(labelsize=10, direction='in', length=5,
                    width=.5, colors='k', grid_color='k', grid_alpha=0.5)

        ax[0].plot(data[:, 0]/1.0e5, data[:, 1], color=color["royal_blue"], linestyle="-", linewidth = lw)
        ax[0].set_xlabel("$r~$[ km]")
        ax[0].set_ylabel("$M/M_{\\odot}$")
        ax[0].set_title("Mass-Distribution")

        ax[1].tick_params(labelsize=10, direction='in', length=5,
                        width=.5, colors='k', grid_color='k', grid_alpha=0.5)
                        
        ax[1].plot(data[:, 0]/1.0e5, data[:, 2], color=color["salmon"], linestyle="-", linewidth = lw)
        ax[1].set_xlabel("$r~$[km]")
        ax[1].set_ylabel("$P/\\epsilon_0$")
        ax[1].set_title("Pressure-Distribution")
        plt.savefig("fig.png")
        plt.show()

#-------------------------------------------------------------------------------------------------

    def single_plot(self, filename, title="R${_{NS}}$ Vs M/M${_{\\odot}}$"):
            try:
                import matplotlib.pyplot as plt
            except:
                print("matplotlib.pyplot ImportError: @ Plot.single_plot")

            try:
                import numpy as np
            except:
                print("Numpy ImportError: @ Plot.single_plot")

            #-------------------------------------------------------------------------------------------------

            plt.rc('text', usetex=True)
            plt.style.use('seaborn-white')

            #-------------------------------------------------------------------------------------------------

            lw = 1.0
            color = {
                "indigo": "#4B0082",
                "royal_blue": "#4169E1",
                "crimson": "#B22222",
                "salmon": "#FA8072",
            }

            #-------------------------------------------------------------------------------------------------

            data = np.genfromtxt(filename)

            _, ax = plt.subplots(1, 1)
            ax.tick_params(labelsize=10, direction='in', length=5,
                            width=.5, colors='k', grid_color='k', grid_alpha=0.5)

            ax.plot(data[:, 0]/1.0e5, data[:, 1], color=color["royal_blue"],
                    linestyle="-", zorder = 1, linewidth=lw)
                    
            ax.scatter(data[:, 0]/1.0e5, data[:, 1], color=color["salmon"],
                     marker='.', zorder = 2, linewidth=lw)

            ax.set_xlabel("R [ km]")
            ax.set_ylabel("M/M${_{\\odot}}$")
            ax.set_title(title)

            plt.show()

#-------------------------------------------------------------------------------------------------
