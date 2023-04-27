# %%
import pandas as pd

df = pd.read_csv("/Users/haoyu/Documents/GitHub/NormalRegPanelMixture/results/powerTest/power.simulate.csv")

# %%
df.columns = ['V0', 'power', 'T', 'N', 'mu', 'alpha', 'sigma']
df = df[['power', 'T', 'N', 'mu', 'alpha', 'sigma']]

# %%
df.pivot(columns=["alpha", "N", "T"], index=["mu", "sigma"], values="power").to_excel("/Users/haoyu/Documents/GitHub/NormalRegPanelMixture/results/powerTest/powerM2.xlsx")

# %%
