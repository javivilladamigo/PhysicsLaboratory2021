# Cloning a repository
To grab a complete copy of another user's repository, use git clone like this:

`$ git clone https://github.com/javivilladamigo/PhysicsLaboratory2021.git`

`# Clones a repository to your computer`

# Fetching changes from a remote repository
Use `git fetch` to retrieve new work done by other people. Fetching from a repository grabs all the new remote-tracking branches and tags without merging those changes into your own branches.

If you already have a local repository with a remote URL set up for the desired project, you can grab all the new information by using git fetch *remotename* in the terminal:

`$ git remote add upstream https://github.com/javivilladamigo/PhysicsLaboratory2021.git`

`$ git fetch upstream`

`# Fetches updates made to a remote repository`


# Merging changes into your local branch

Merging combines your local changes with changes made by others.

Typically, you'd merge a remote-tracking branch (i.e., a branch fetched from a remote repository) with your local branch:


`$ git merge upstream/main`

`# Merges updates made online with your local work`

# Pulling changes from a remote repository

`git pull` is a convenient shortcut for completing both git fetch and git merge in the same command:

`$ git pull upstream main`

`# Grabs online updates and merges them with your local work`

Because `pull` performs a merge on the retrieved changes, you should ensure that your local work is committed before running the `pull` command. If you run into [a merge conflict](/github/collaborating-with-pull-requests/addressing-merge-conflicts/resolving-a-merge-conflict-using-the-command-line you cannot resolve, or if you decide to quit the merge, you can use `git merge --abort to take the branch back to where it was in before you pulled.
