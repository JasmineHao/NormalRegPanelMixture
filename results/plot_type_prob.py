# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import itertools
from textwrap import wrap

# %%

models = ["k_dummy", "k_dummy_imex", "k_dummy_ciiu", "k_dummy_imex_ciiu"]

models = ["k_dummy_imex", "k_dummy_imex_ciiu"]

industries = ['Food products','Fabricated metal products','Textiles']
industry = industries[0] 
model = models[0]

# %%
for model in models:
    for industry in industries:
        
        data_plot = pd.read_csv( os.path.join(r'C:\Users\haoja\Documents\GitHub\NormalRegPanelMixture\results\Empirical\Additional', 
                                model + "_" + industry + "_type_prob.csv"))

        # data_plot = pd.read_csv( os.path.join(r'C:\Users\haoja\Documents\GitHub\NormalRegPanelMixture\results\Empirical\Additional', 
        #                         model + "_" + industry + "_statistics.csv"))
        
        print(model + "_" + industry + "_type_prob.csv")
        data_plot = data_plot.sort_values("group_num")
        data_plot = data_plot[[each for each in data_plot.columns if 'comp' in each]]
        data_plot.columns = [int(each.strip('comp.')) for each in data_plot.columns]

        # Define the categories
        capital = ["low capital", "high capital"]
        importer_status = ["non-importer", "importer"]
        exporter_status = ["non-exporter", "exporter"]

        # Generate all combinations using itertools.product
        combinations = list(itertools.product(capital, importer_status, exporter_status))
        # combinations = [",".join(each) for each in combinations]

        report_val = data_plot.T

        if report_val.shape[1] == 8:
            report_val.columns = combinations
        else:
            report_val.columns = capital
        

        # Example data: 9 types, each with 8 sub-types
        types = np.arange(len(report_val)) + 1
        subtypes = np.arange(len(report_val.columns))
        
        values = report_val  # Random values for illustration
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(10, 6))

        # Define the bar width and x locations for each sub-type
        bar_width = 0.1
        x_positions = np.arange(len(subtypes))

        # Plot for each type
        for i in range(len(types)):
            # Calculate positions for the current type's subtypes
            x_offset = i * bar_width
            ax.bar(x_positions + x_offset, values.loc[i+1,], width=bar_width, label=f'Type {types[i]}', capsize=5)

        # Customize plot
        ax.set_xlabel('Sub-types')
        ax.set_ylabel('Values')
        ax.set_title('Observed Types and Unobserved Types: ' + industry  + ' ' + model )
        ax.set_xticks(x_positions + (bar_width * (len(types) - 1)) / 2)  # Set the x-tick labels
        if len(values.columns) == 8:
            ax.set_xticklabels(['\n'.join(each) for each in values.columns])
        else:
            ax.set_xticklabels(values.columns)
        ax.legend()

        # Show plot
        plt.tight_layout()
        plt.show()

        
        fig_path = r'C:\Users\haoja\Documents\GitHub\NormalRegPanelMixture\results\Empirical\figure'
        fig.savefig(os.path.join(fig_path, model + '_'  + '_'.join(industry.lower().split(' '))  + '_prob.png'))


# %%
