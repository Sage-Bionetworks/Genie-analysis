
## Contributing

We welcome all contributions! This repository is structured such that each analysis or project will have its own folder. Please head to [issues](https://github.com/Sage-Bionetworks/Genie-analysis/issues) to either file any analysis requests/errors or find a task you want to assist with.  Make sure to assign yourself the task if you decide to work on it.


### Fork and clone this repository

See the [Github docs](https://help.github.com/articles/fork-a-repo/) for how to make a copy (a fork) of a repository to your own Github account.

Then, [clone the repository](https://help.github.com/articles/cloning-a-repository/) to your local machine so you can begin making changes.

```
git clone https://github.com/<your-gh-user>/Genie-analysis.git
```

Add this repository as an [upstream remote](https://help.github.com/en/articles/configuring-a-remote-for-a-fork) on your local git repository so that you are able to fetch the latest commits.

```
cd Genie-analysis
git remote add upstream https://github.com/Sage-Bionetworks/Genie-analysis
```

On your local machine make sure you have the latest version of the `main` branch:

```
git checkout main
git pull upstream main
```

### Commit and Create a Pull Request

Create your own branch from the `main` branch.  Create your own analysis folder or choose to put your analysis code in an existing correct project or analysis folder. Here I am naming my branch `analysis-a`, but it can be something more descriptive:

```
git checkout -b analysis-a
git push --set-upstream origin analysis-a
```

[Commit the scripts or changes](https://guides.github.com/introduction/git-handbook/#basic-git) you would like and create a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) against this respository.  Your contribution will be reviewed and merged.
