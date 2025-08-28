import re
import pandas as pd

# Regex para capturar Iter, N.size, A.size
pat = re.compile(
    r'Iter T=(?P<Iter>\d+) \| '
    r'N\.size\(\)=(?P<Nsize>\d+) .* '
    r'A\.size\(\)=(?P<Asize>\d+)'
)

data = []
with open('C:/Users/Giova/OneDrive/Escritorio/Mining/alg/LG/debug.log','r') as f:
    for L in f:
        m = pat.search(L)
        if m:
            data.append(m.groupdict())

df = pd.DataFrame(data).astype(int)
print(df.head())
# Opcional: guardar a CSV
df.to_csv('debug_summary.csv', index=False)
