---
layout: post
title: My New Notebook
date: '2018-11-08'
tags:
  - computing
  - wordpress
  - jekyll
  - onsnetwork
categories:
  -
---
Here we are - my new notebook!

I finally gave up on [my onsnetwork notebook](onsnetwork.org/kubu4). I hit a few rough patches recently where the site was down and it was enough to push me over the edge. I needed a simpler way of having an offline version of my notebook available so that I could view the contents quickly and easily even when there are network outages.

So, I moved things over to here: GitHub Pages. Everyone else in the lab has been using these for quite some time.

It took me a bit of time to get everything configured, as I wanted to use a theme, instead of the default "jekyll-now" theme that everyone has has used. This theme ([Basically Basic](https://mmistakes.github.io/jekyll-theme-basically-basic/)) had all of the features that I felt I needed:

- Built-in search

- Tag organization

- Project (i.e. Categories) organization

- Display post date alongside the post excerpt

Getting things up and running just to begin posting wasn't terribly difficult. I set things up using the "remote themes guide" found in the [Basically Basic](https://mmistakes.github.io/jekyll-theme-basically-basic/) documentation.

However, trying to get all my posts moved over from [my onsnetwork notebook](onsnetwork.org/kubu4) was a bit more challenging. I exported my notebook (which is Wordpress) to an XML file and used the [Exitwp](https://github.com/thomasf/exitwp) program to convert to Jekyll (the backend software that runs on GitHub Pages to generate the site) markdown files. The conversion hit some snags that I still haven't resolved:

![](https://raw.githubusercontent.com/RobertsLab/sams-notebook/master/images/screencaps/20181107_001.png)

I'll have to manually transfer the content of those two posts.

I also needed to convert URLS within posts that linked to other notebook posts to link to other Jekyll posts and not back to [my now-obsolete notebook](onsnetwork.org/kubu4).

I wrote the following script to do that conversion:

- [jekyll_url_converter.sh (GitHub)(https://github.com/kubu4/Scripts/blob/master/bash/jekyll_url_converter.sh)

Although that script looks short and sweet, it took an inordinate amount of time to make it functional (I really got hung up on the regex in the ```grep``` command). I had to post to [Stack Overflow](https://stackoverflow.com/questions/53213379/grep-o-multiple-occurrences-of-variable-string-in-same-line) to get it figured out.
