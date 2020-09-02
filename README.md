# icesat2-book

### Project description
Code for building the ICESat-2 Jupyter Book.<br><br> Link to jupyter book: (https://nicolejkeeney.github.io/icesat2-book/home.html)[https://nicolejkeeney.github.io/icesat2-book/home.html]

### Steps to build book 
 1) Delete **build** and **pycache** folders
 2) Update github repo with any changes 
 3) **cd** out of local book directory 
 4) In terminal, run: **jb build icesat2-book** 
 5) **cd** into local book directory 
 6) In terminal, run: **ghp-import -n -p -f _build/html** 
