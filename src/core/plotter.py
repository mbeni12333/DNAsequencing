import matplotlib.pyplot as plt

class plotter:
    def __init__(self, xsize=150, ysize=150):
        #size, figure
        return
    def plot(self, x, y):

        plt.plot(x, y)

    def saveplot(self, filename=""):
        return
    def show(self):
        plt.show()

if __name__ == "__main__":
    p = plotter()
    x = list(range(10))
    y = [i*i for i in x]
    y2 = [i**(1/2) for i in x]
    p.plot(x, y)
    p.plot(x, y2)
    p.show()
