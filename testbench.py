import tkinter as tk
from analysis_expt import xpm_analysis



root = tk.Tk()
root.minsize( width=760, height=530)
ws = xpm_analysis(root)
root.mainloop()
