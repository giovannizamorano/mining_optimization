import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

df = pd.read_csv("C:/Users/Giova/OneDrive/Escritorio/Mining/alg/LG/pit_blocks.csv")

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df['X'], df['Y'], df['Z'], c=df['Value'], cmap='viridis', s=10)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Pit Final (Lerchs-Grossmann 3D)')

plt.show()