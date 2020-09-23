<h1 align="center"> Rules of Operation </h1>
<p align="center"> <em>For a Clean, Efficient, and Transparent Collaborative Workflow</em></p>

> This file is intended to serve as a set of guidelines for us to follow that should help reduce confusion, enforce consistency, and facilitate the development process. It is based on B. Staton's past experience, the experience of the modeling team in a practice collaborative project, and various websites that document best practices for Git/GitHub.
>
> **By pushing commits to this repository, you are agreeing to follow these rules to the best of your ability**. Of course, mistakes happen and we promise to be understanding, but please strive to follow these rules to the greatest extent possible.
>
> Suggested changes may be made with a `doc-rule-suggestions` branch and pull request. These guidelines may be updated as we move forward.
>
> _PLEASE READ IT CAREFULLY PRIOR TO MAKING ANY COMMITS_

## Git/GitHub-Related Topics

**NOTE**: It is recommended to use the [GitHub Desktop](https://desktop.github.com/) app for interacting with the Git/GitHub aspects of your work on this project. See [here](https://docs.github.com/en/desktop) for user documentation.

### Workflow

We are using the [**Feature Branch Workflow**](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow). This means:

* There is **one remote repository** - no forks.
* Other modelers will be added as **collaborators and they clone** the remote to a local copy.
* The **master branch** is viewed as the most current "deploy-able" branch - it represents the best available version of the model at the time.
* All work is done in **feature branches** off of master.
  * **NEVER** commit/push directly to master without first discussing it with the repo owner. This will prevent merge conflicts.
* **Pull requests** are submitted and reviewed prior to merging them into master.
  * As the repository owner, it is B. Staton's responsibility to be the primary reviewer and merger.

### General

* **Pull** from the remote frequently.
  * Pull before creating a branch.
  * This will help prevent merge conflicts.
* **GitHub issue tracker** should be used to document, discuss, and delegate main tasks.
* **Markdown syntax** should be used in all correspondence on GitHub. E.g., if you are referring to code, make it `monospaced`.
  * An overview of GitHub-flavored markdown can be found in the GitHub guide [Mastering Markdown](<https://guides.github.com/features/mastering-markdown/>).
* **Test your code** before committing, and especially before pushing.

### Commits

* Commits should be **atomic** - _one change per commit_.
  * This is very important for code review.
  * Also rolling back changes (less familiarity here)
  * Excellent article [here](https://www.freshconsulting.com/atomic-commits/) for more details on why this is important and some basic rules.
  * This **does not mean** that each line of code you change should be a separate commit.
  * Rather, it **does mean** that each batch of changes that accomplish a single task should be a commit.
* Examples of atomic commits (example commit message summary lines given):
  * Add data files X, Y, Z
  * Add data preparation for juvenile survival
  * Rename variable `abc` to `ABC` throughout
  * Edit code for figure
  * Add function to perform X
* You do not need to commit all changes at once.
  * GitHub Desktop is great for selecting which lines/files you wish to commit.
* Don't commit just to get a clean branch to allow switching branches.
  * Use the **[stash](https://github.blog/2019-06-05-github-desktop-expands-to-support-stashing-and-rebasing/?utm_campaign=1559705923&utm_medium=social&utm_source=linkedin&utm_content=1559705923#stashing)** feature instead.
* Commits should be accompanied by a **good commit message**.
  * Often the default is inadequate ("Update FileXYZ"); unless it is a README.
  * [Here](https://chris.beams.io/posts/git-commit/) is an article about the importance of and suggestions for writing good commit messages - no need to worry about _exact_ details here.
  * The **message summary** field should contain a title that states the task accomplished by the commit. It should typically start with a verb ("Add", "Update", "Fix", "Refactor", "Delete", "Edit", etc.) - generally 50 characters or less.
  * The **message description** field elaborates on what was changed and why. Bullet points with "*" or "-" are suggested.
  * If a commit addresses a particular issue, you should reference it using `#`.
  * [Here](https://imgs.xkcd.com/comics/git_commit_2x.png) is an entertaining comic representing what we are **trying to avoid**.

### Branches
* Recall all changes should be made on branches off master.
  * They are merged via **pull request**.
* Branches should be relatively small and constitute one feature at a time.
* Branches may (**should**) be deleted after they are merged.
  * **NOTE**: Branches deleted from the remote will remain in local repos until you manually delete them.
* Branches should start with a prefix denoting what that branch does. Examples:
  * `feature-add-bio-data`
  * `bugfix-fix-spelling`
  * `docs-add-meeting-notes`
  * `refactor-switch-variable-names`
* In general, there should be one collaborator working on one branch.
  * Two or more people pushing commits to a branch can create conflicts in the history.
  * Frequent pulling might allow this to work without issues.
  * If it is necessary for two or more people to push commits to a single branch, then **close communication** should be had to ensure their local repos are in sync.
* If master moves forward from the point where your branch was created, your branch should generally be rebased prior to merging it. 
  * This keeps the commit history clean.
  * This may require some experimentation.

## Code-Related Topics

* Spaces should not be used in file names and all words should be lowercase, except the extension (e.g., `.R`).
  * Use `-` **not** `_` to separate words in a file name.
  * E.g., `jags-model-code.R` not `jags_model_code.R`
* Object names should be concise, but descriptive.
    * Do not use things like `a`, `b`, etc. to represent important objects.
    * `_` is preferred to `.` in object names, e.g., `juv_surv`, not `juv.surv`
* All R packages used by **any code** in this project should be placed in the file `00-packages.R` using `library(pkgName)`.
  * This script should be sourced at the top of R files using `source("00-packages.R")`.
  
  * If you include a new package into this file, include a brief statement about what it is used for and how to obtain it, e.g.,
    ```R
    library(postpack) # summarizing posterior samples; install.packages("postpack")
    library(jagsUI)   # calling JAGS from R; install.packages("jagsUI")
    ```
  * To the extent practicable, we should strive to minimize dependencies
    * The more packages we depend on, the more ways our code can break due to factors outside of our control.
* All functions should be placed in the `01-functions` directory, and placed in the appropriate file
  * `data-fns.R`: functions used in preparing raw data files for use in the rest of the analysis
  * `jags-fns.R`: functions used in calling JAGS - these may include functions for setting initial values, calling JAGS, bundling data, etc.
  * `plot-fns.R`: functions for plotting
  * `sim-fns.R`: functions for simulating populations
  * `util-fns.R`: basic utility functions that serve some menial task
  * If your function doesn't belong in one of these categories, create a new file in `01-functions` using the same syntax as above. 
  * We can load these functions into any R session using (this will go at the top of most scripts, along with `source("00-packages.R")`)
     ```R
     # load all functions into session
     invisible(sapply(list.files("01-functions", full.names = T), source))
     ```
* All file paths should be **relative to the project directory**. This will ensure the code is reproducible regardless of who is running it on their local machine. When you work on code for this project, make sure the RStudio project is active in your R session.
* All code should be **responsibly annotated with comments**. No need to be excessive (e.g., every line), but it should always be clear what purpose all portions of the code serve.
* Any files you do not wish to commit should be placed in your personal `98-scratch` directory. **All files in this directory are ignored by Git**. This means its contents are unaffected by branching and will never be committed. You will need to create this directory yourself after cloning the remote repo. 
* Some prefer to use `=` for assignment in R rather than `<-`, but given these serve an equal purpose, either may be used.
  * When a "final version" of the repo is ready, we may perform a replace for consistency throughout
* Use limited (if any) reliance on the "tidyverse"
  * Exceptions for data exploration with ggplot2
* Subsets should never be hardcoded with element numbers
  * Use column names, logicals, etc.
  * E.g., to extract spawner data for a particular population:
    ```R
    dat[dat$pop == "the_pop","spawners"]  # YES
    dat[1:20,10]                          # NO
    ```
