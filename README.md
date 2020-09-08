# icesat2-book

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nicolejkeeney/icesat2-book/master)

### Project description
Code for building the ICESat-2 Jupyter Book.<br><br> Link to jupyter book: [https://nicolejkeeney.github.io/icesat2-book/](https://nicolejkeeney.github.io/icesat2-book/)<br><br>This repository is configured to run in Binder.



### Activating the environment 
Run **conda env create -f environment.yml** in terminal.<br>To activate the environment, run **conda activate is2-book-env**



### Steps to build book 
 1) Update github repo with any changes
 2) Activate virtual notebook
 3) **cd** out of local book directory 
 4) In terminal, run: **jb build icesat2-book** 
 5) **cd** into local book directory 
 6) In terminal, run: **ghp-import -n -p -f _build/html** 
 
