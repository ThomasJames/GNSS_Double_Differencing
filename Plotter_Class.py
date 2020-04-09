import seaborn as sns
import matplotlib.pyplot as plt

class HeatMap:

    def __init__(self, matrix, title):
        self.matrix = matrix
        self.title = title
        self.file_type = ".png"
        self.colour = 'Pastel1'

    def output_png(self):
        sns.heatmap(self.matrix,
                    annot=True,
                    cbar=False,
                    xticklabels=False,
                    yticklabels=False,
                    cmap=self.colour)
        plt.title(str(self.title))
        plt.savefig("Matrix_Output/" + str(self.title) + self.file_type)

    def matrix_plotter(self):
        sns.heatmap(self.matrix,
                    annot=True,
                    cbar=False,
                    xticklabels=False,
                    yticklabels=False,
                    cmap=self.colour)
        plt.title(str(self.title))
        plt.show()


