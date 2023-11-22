# Balanced: removal of 10% from top and bottom.

setting <- "balanced"
ycut_u = 0.9
ycut_l = 0.1
source('Code/imputation_sim_simplified_MNAR.R')

# Top: removal of 10% from top

setting <- "top"
ycut_u = 0.9
ycut_l = 0 
source('Code/imputation_sim_simplified_MNAR.R')

# Bottom: removal of 10% from bottom

setting <- "bottom"
ycut_u = 1
ycut_l = 0.1
source('Code/imputation_sim_simplified_MNAR.R')

# Unbalanced: removal of 10% from top and 5% from bottom

setting <- "unbalanced"
ycut_u = 0.9
ycut_l = 0.05
source('Code/imputation_sim_simplified_MNAR.R')