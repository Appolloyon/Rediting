Rediting
========

This repository contains several Python scripts to calculate various metrics of
transcript editing in biological data. These include characterizing individual
editing events by their ability to modify sequence properties, comparison of
sequences to established references to determine biological trends of editing,
and running basic simulations to assess editing properties.

All scripts in this repo are authored by Christen M. Klinger
email: cklinger@ualberta.ca

Setup
=====

These scripts should run on any computer with a working Python installation,
though they have not been explicitly tested on Linux or Windows platforms.

If you do not have a working Python installation, obtain this first, e.g. by
following the recommendations at:

https://wiki.python.org/moin/BeginnersGuide/Download

Once this is done, you can either download the code from this page directly,
or (after installing git), simply for this repo.

Alternatively, create a virtualenv with python version ~2.7.6 called rediting

```
    $ virtualenv ~/envs/rediting --python python2.7 
    $ source ~/envs/rediting/bin/activate
```

Assuming you have pip installed, install the requirements:

```
    $ pip install -r requirements.txt --upgrade
```

The scripts can then be run on your system, either by providing the full
relative path to the directory holding the code on your system, or by adding
this directory to your .bash_profile (at least on Unix systems).

Running
=======

All scripts in this folder have some built-in documentation. To view this,
simply call the script using either the -h or --help command line flags.

Some things to keep in mind:
- All genomic and reference sequences should be in the proper reading frame
- The scripts require pre-aligned sequences. This is a slight burden to the
user, but also allows the use of many alternatives (e.g. MAFFT, MUSCLE, etc.)

A Note to Users
===============

These scripts are provided as-is, in the hopes that they may be useful. However,
they come with no warranty of any kind, explicit, implicit, or otherwise. If you
wish to use these scripts, and run into issues, please do some initial research
to determine whether it is an easily fixable issue before emailing the author (I
will do my best to respond to any inquiries, but I can make no guarantees).

A Note to Developers
====================

These scripts were coded while I was still learning various features, not only of
Python, but of development in general. Hence, the code is still messy and the
implementation is not ideal. I have plans to revisit this in future, time permitting,
but in the meantime further development/extension of the code is encouraged by
any interested parties.
