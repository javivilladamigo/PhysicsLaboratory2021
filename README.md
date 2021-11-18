# PhysicsLaboratory2021

## Git instructions

### Setting up a local repository

The following instructions need to be followed any time a new local repository is created. If you are working in a location where such repo already exists, what follows doesn't need to be repeated every time.

   * Clone this repository (i.e. create a local repository cloned from this remote repository)

   `git clone https://github.com/javivilladamigo/PhysicsLaboratory2021.git`

   A new directory will appear in your current directory. Get into it:

   `cd PhysicsLaboratory2021/`

   * Define the central PhysicsLaboratory2021 repo as the upstream repository:

   `git remote add upstream https://github.com/javivilladamigo/PhysicsLaboratory.git`

   * Check that the previous commands succeeded:

   `git remote -v`

   * Get (fetch) the updates that have been done on the remote repository:

   `git fetch upstream`

  * The default branch is `main`. This is the branch you should point the pull request to. You should now create your development branch where you can edit and run the code. In order to set up a proper development cycle, you must create a branch (in the example below called `dev`) that *tracks* `upstream/dev`:

   `git branch -vv`

   `git checkout -b dev upstream/dev`

### Standard development cycle

   * Before starting with the development you could check whether the upstream repository has been updated with respect to your forked version (that's likely to be the case prior to every lab class). If it had, then merge the changes into your main:

   `git checkout main`
   
   `git fetch upstream`

   `git merge upstream/main`
   
   * Then synch your development branch (especially in the case your pull request has been recently approved):
   
   `git checkout dev`

   `git fetch upstream dev`

   `git merge upstream/dev`

   Be careful that the git syntax is inconsistent between fetch and merge. In the former you should use the whitespace to separate the repository and the branch name, in the latter you should use the slash character.

   * The idea is that your main always reflects `upstream/main`, i.e. it keeps a local copy of the reference code as a starting point for your developments (i.e. solving the assigned problems). Note that in order to update your repository on GitHub, you need to push the local version to your remote repository.

   * You may also need to get the updates from the main, i.e. need to merge the main:

   `git merge main`

   * Before starting to edit on the machine that you are using, type the follow command in order to update the directory with the last changes:
  
   `git pull`

   * Now develop some code. Image you create a `<NewFile>`. Add the file to your local repository and stage it for commit (to unstage a file, use `git reset HEAD <NewFile>`)

   `git add <NewFile>`

   * Commits the (tracked) changes you made to the file(s) and commit them local repository on github

   `git commit -m "Add existing file"`

   (what follows after `-m` is a comment to keep track of the reason of the commit)

   * Now propagate (push) your local changes to your remote repository on github (`origin`)

   `git push origin dev`

   * When appropriate, propagate your development also to the main repository (upstream). For that you need to go for a pull request, which is done from GitHub. Pay attention to set the correct starting and destination branches. This should be done from your dev branch, pointing towards the main repository. After that I will accept the pull request and commit the changes.
