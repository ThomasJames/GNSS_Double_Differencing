import matplotlib.pyplot as plt
from statistics import mean
from scipy.stats import zscore
import seaborn as sns


x = [4929605.54130538,  4929605.5410819, 4929605.53911663]
y = [-29123.82730074, -29123.82767789, -29123.82739046]
z = [4033603.93206799, 4033603.93099312, 4033603.93271809]


# Get the mean of the 3 epochs.
x = mean(x)
y = mean(y)
z = mean(z)


# collect the residuals from the three epochs
b = [9.29757378e-01, 2.52782310e-04, 9.89184036e-01, 2.90583952e-01, 9.64209316e-01, -1.13994574e-01,7.92398563e-02,
     0.94013581, 0.00273867, 0.97883069, 0.28919151, 0.97179478, -0.10708197, 0.08768272, 0.93938436, 0.01976735,
     0.98053831, 0.29998728, 0.97263572, -0.10390881, 0.06579685]

# Get the Z-Scores of the b vector
b_z = zscore(b)

# Plot the Z-Scores
sns.distplot(b_z, bins=200, kde=False)
plt.xlabel("Z Scores")
plt.ylabel("Frequency")
plt.title("Z Scores of Observed - Computed")
plt.savefig("Graphs/Z Scores of Observed - Computed") 
plt.show()








