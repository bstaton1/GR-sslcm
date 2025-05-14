# SET WORKING DIRECTORY TO GR-SSLCM.RPROJ LOCATION
proj_dir = this.path::this.proj()

# this directory
this_dir = file.path(proj_dir, "03-post-process/ms-content")

# scenario name
if (!exists("scenario")) scenario = "test_vshort"

# set the proportion of MCMC samples to retain
keep_percent = 1

# set the file extension of figures
if (!exists("fig_type")) fig_type = "png"

# load the model info
source(file.path(this_dir, "load-model-info.R"))

# rebuild all figure and table content from source
my_cat("Creating MS Content:", "No Fall Migrants (Tab 1)", first = TRUE)
source(file.path(this_dir, "make_table_no-fall.R"))

my_cat("Creating MS Content:", "Obs vs. Fit (Fig 3)")
source(file.path(this_dir, "make_figure_obs-v-fit.R"))

my_cat("Creating MS Content:", "Rates Comparisons (Figs 4 & 6)")
source(file.path(this_dir, "make_figure_rate-compare.R"))

my_cat("Creating MS Content:", "WUL vs. Capacity (Fig 5)")
source(file.path(this_dir, "make_figure_capacity-v-wul.R"))

my_cat("Creating MS Content:", "Relationships (Fig 7)")
source(file.path(this_dir, "make_figure_relationships.R"))

my_cat("Creating MS MS Content:", "Rho Comparisons (Fig 8)")
source(file.path(this_dir, "make_figure_rho-compare.R"))

my_cat("Creating MS Content:", "WUL Change Scenarios (Fig 9)")
source(file.path(this_dir, "make_figure_wul-change.R"))
my_cat("Creating MS Content:", "DONE")