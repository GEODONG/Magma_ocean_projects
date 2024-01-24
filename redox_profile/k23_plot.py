import pandas as pd
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

# 1. Load the dataset from the Excel file.
data = pd.read_excel("k23.xlsx",sheet_name="fO2_MO_adiabat")

# 2. Extract the columns
X = data[" Fe3+/ΣFe=0.03"].values.reshape(-1, 1)
y = data["P (GPa)"].values

# 3. Use polynomial regression to fit the data.
degree = 3
poly_features = PolynomialFeatures(degree=degree)
X_poly = poly_features.fit_transform(X)

model = LinearRegression().fit(X_poly, y)

# Predict
X_seq = np.linspace(X.min(), X.max(), 300).reshape(-1, 1)
y_pred = model.predict(poly_features.transform(X_seq))

# 4. Plot the data and the fitted curve
#plt.scatter(X, y, color='blue', label='Actual Data')
plt.plot(X, y, color='blue', label='Actual Data')
#plt.plot(X_seq, y_pred, color='red', label=f'Fitted Polynomial (degree={degree})')
plt.xlabel("Fe3+/ΣFe")
plt.ylabel("P (GPa)")
plt.title("Polynomial fit of Fe3+/ΣFe vs P (GPa)")
plt.legend()

plt.gca().set_ylim(plt.gca().get_ylim()[::-1])

plt.savefig("./k23.png")
#plt.show()
#" Fe3+/ΣFe=0.03"