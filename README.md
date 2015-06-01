# INSTALLATION - SET-UP

In order to use Makepipe to create a new pipeline project, first clone this
repository somewhere on you computer. Then use the `build.sh` script to
generate a tarball that only packup thinks you need to deploy you project:

    ./build.sh # This will create a Makepipe-X.XXX.tar.gz tarball

Then all you need to do is unpack the tarball on your project directory and
start writing the YAML pipeline config.

    cd my_project
    cp /XXX/XXX/Makepipe-X.XXX.tar.gz .
    tar -xzf Makepie-X.XXX.tar.gz
    vi pipeline.yml

# USER DOCUMENTATION

Makepipe helps the creation of simple to sofisticated pipeline relying on
Makefile technology without getting your hands dirty. It has been concepted to
create collections of re-usable bricks that can be shared among multiple
projects. It has very low dependancies (Perl core and make), and it is very
quick to set-up.

Create a pipeline with Makepipe is simple as hell, you just need to create a
simple YAML description of the pipeline and we do the rest, that is generating
a fully fonctionnal Makefile version of the pipeline, ready to be "maked".

First let see how to create this YAML pipeline.

# WRITING YAML PIPELINES

YAML is a syntax language used to describe data's in a human friendly way
(unlike XML that is hard to write, impossible to read, and not that easy to
parse. 

In order to create your first pipeline, start by creating and empty YAML file,
with an `.yml` extension.

## Collections

A YAML document is made of collections, lets define our first collection that
will hold our samples names

    ---
    samples:
      - name: test1
      - name: test2

We can also add another collection that will hold some configuration variables

    ---
    samples:
      - name: test1
      - name: test2
    config:
      - genome: GRCh38
      - reads_dir: raw_reads

You can define as many collections as you want, and gives them the name you
desire.

## Brick instances

### My first brick instance

You can also define some special collection that we call "brick instances".
These special collection will be interpreted by Makepipe in order to created
process that will be exectued in the pipeline.

    my_first_brick: 
      brick: GENERIC 
      config: 
        COMMAND_LINE: "ls"

This first brick install creates a process called "my first brick", wich is an
instance of the brick GENERIC (located in the bricks/ folder). This instance of
GENERIC declare is configured to execute the command line `ls`. If you generate
a make file from this pipeline, you'll see that the process "my first brick"
will be executed and print the results of `ls` on the standard output

    > ./makepipe build my_pipeline.yml > Makefile
    > make my_first_brick
    make[1]: Entering directory '/home/jaudoux/Dev/makepipe'
    ls
    bricks    Makefile  test.yml   makepipe
    make[1]: Leaving directory '/home/jaudoux/Dev/makepipe'make

### Run only a given brick

Makepipe has automatically place your block instance under the "all" rule, so
it executed automatically when make is called without arguments. You can also
specifically execute this block instance by running `make my_first_brick`.

### Makepipe auto-regeneration

When the Makefile if first generated, makepipe places a rule that will check if
your YAML pipeline has been modified since the Makefile was generated in order
te re-generate the Makefile before running it. In other word, do not care about
running the `makepipe build`, once the Makefile is constructed the only command
you need to do run is `make`.

### My first pipeline

If we look at the man page of the GENERIC brick by typing `./makepipe brick
GENERIC`, we can see that there is two other parameter that we can give to
GENERIC : `INPUT_FILE` and `OUTPUT_FILE`.

We can now create our first pipeline with two bricks that will be chained by IO
files.


    brick1:
      brick: GENERIC
      config:
        OUTPUT_FILE: test.txt
        COMMAND_LINE: "ls > test.txt"
    brick2:
      brick: GENERIC
      config:
        INPUT_FILE: test.txt
        COMMAND_LINE: "cat test.txt"

The link between these two instances is only described by the file test.txt
that is the output of the first brick and the input of the second one. You do
not need to specify anything else, GNU make will do the rest.

The first time the pipeline will be run, the brick1 will be called, creating
the test.txt file needed by brick2, then brick2 will cat the file to the
standard output. If you call the pipeline a second time, only brick2 will be
called this test.txt already exist.

## BRICKS COLLECTION

Until there, we have only declared instances of the GENERIC brick, but there is
many others and you can even create a re-usable brick that you want.

Basically a brick is a pre-formated Makefile that delares some variables and
define some rules.  When you instatiate a brick, the Makefile is read and
variables are replaced, and integrated in the output Makefile. Each brick
define its own variables, and you must read the brick documentation before
creating an instance.

The subcommand `brick` will list all available bricks that are located under
the `bricks/` directory.

    ./makepipe brick
    Availables bricks are:
      CHIMCT
      CORE
      CRAC
      CRACTOOLS_EXTRACT
      FEATURECOUNTS
      GENERIC
      GZIP

Then you can call the subcommand `brick` with the brick name as argument to get
information about it. For example `./makepipe brick CRAC`.

You'll see that some variables are automatically computed, some others have
default values, but some might be required to create the brick instance. Those
variables are assigned to the value "undef", meaning that you need to give a
value for this variable under the "config" entry of the brick instance.

A brick can be instiate multiple times like we have done with GENERIC brick in
the previous example.

## Use YAML collections and variables {{..}}

In the first part we have learn to declare collections that can be used to
store all kind of informations (including bricks configuration), but how can we
use/access these informations?

In Makepipe you can use a syntax inspired by the broadly known ["moustache
syntax"](https://mustache.github.io), by placing the variable name you want to
acces inside double braces tag.

Because an example is always easier to understand than a bunch of litterall
explainations, here's one:

    directory: /usr/bin
    list_bin:
      brick: GENERIC
      config:
        COMMAND_LINE: "ls {{directory}}"

In this example the brick instance (called `list_bin`) of GENERIC have specify
a command line that list the content of the directory declared at the begining.
This can be helpfull when different bricks share some similar variables.

Remember to place all expressions that uses curly braces into quotes, otherwise
it will be interpreted by YAML parser as a hash structure!

The syntax can be extended to access variables that are stored in hash
collections by using a dot as separator :

    config:
      directory: /usr/bin
      ls_command: ls
    list_bin:
      brick: GENERIC
      config:
        COMMAND_LINE: "{{config.ls_command}} {{config.directory}}"

Will see later how to access array collections with some kind of loop structures.

## Loops

Sometimes, you want to repeat a given process on multiples values, making some
kind of loops over a brick instance. This is allowed in Makepipe by the `loop`
keyword, that enable multiple iterations of a brick instance on an array
collection.

For example, if we want to list the content of multiple folders we can do the
following brick delaration:

    folders:
      - /usr
      - /tmp
      - /home
    list_folders:
      brick: GENERIC
      config:
        COMMAND_LINE: "ls {{item}}"
      loop:
        list: folders 

In this case the brick `list_folder` will be intantiate multiple times by loop
over the list "folders". At each iteration a special collection "item" can be
used to get the current list item value.

We can also use more complicated array collections, were each item is a hash :

    folders:
      - name: "temporary"
        url: /tmp
      - name: "home"
        url: /home
    list_folders:
      brick: GENERIC
      config:
        COMMAND_LINE: "ls {{item.url}}"
      loop:
        list: folders 
        id: name

In this example we can access the value of the named entry "url" of each item.
We also set the loop id to "name" in order to have nice target names in the generated Makefile.
This allows us to run the `list_folder` brick for only one of the list item :

    > make list_folders_home

This make call will only exectute the block for the /home directory.

## Collection functions

The following functions can be applied on a collection :

### extract()

The extract function is able to extract the value of a given variabe in a collection
of items. For example, we can use the previous example and provide a single rule
that will list all the folders wil a single `ls` call :

    folders:
      - name: "temporary"
        url: /tmp
      - name: "home"
        url: /home
    list_folders:
      brick: GENERIC
      config:
        COMMAND_LINE: "ls {{folders.extract('url')}}"

The Makefile generated from this YAML pipline will exectute the single command
`ls /tmp /home`.

## Exporting

Sometimes you want to save some variables that have been declared in a brick
and store it in a collection. This is allowed by the "export" keyword.

For example imagine that we re-use our first pipeline but we modify it to make
it run on multiple directories:

    folders:
      - name: "temporary"
        url: "/tmp"
      - name: "home"
        url: "/home"
    # Run ls command on a list of folders and place the results in a file
    brick1:
      brick: GENERIC
      config:
        OUTPUT_FILE: "{{item.name}}_ls_file.txt"
        COMMAND_LINE: "ls {{item.url}} > {{this.OUTPUT_FILE}}"
      loop:
        list: folders
        id: name
      export:
        - value: "{{this.OUTPUT_FILE}}"
          to: item.ls_file
    # Print files created by brick1
    brick2:
      brick: GENERIC
      config:
        INPUT_FILE: "{{item.ls_file}}"
        COMMAND_LINE: "cat {{item.ls_file}}"
      loop:
        list: folders
        id: name

The brick1 is now looping over the folders and the file where the ls command is
outputed is stored into each item under the `ls_file` name. The brick2 also
loop over the "folders" collection and uses the `ls_file` new entry created by
brick1.

Remind that you can export as many values as you want.

## EXAMPLE: RNA-Seq pipeline - WIP

In this example we will see how we can pull some bricks together in order to
create a fully fonctionnal RNA-Seq pipeline with makepipe, CRAC and
FeatureCounts.

    ---
    config:
      nb_threads: 4
      nb_jobs: 2
      genome_name: GRCh38

    samples:
      - name: test1 
        type: TOTO
      - name: test2 
        type: TATA

    crac:
      brick: CRAC
      #execution: serial
      config:
        INDEX: "/data/indexes/crac/{{config.genome_name}}"
        OPTIONS: --no-ambiguity
        KMER_LENGTH: 22
        READS: "raw/{{item.name}}_1.fastq raw/{{item.name}}_2.fastq"
        OUTPUT_PREFIX: "mapping/crac/{{item.name}}"
        NB_THREADS: "{{config.nb_threads}}"
      loop: 
        list: samples
        id: name
      export: 
        - value: "{{this.OUTPUT_FILE}}"
          to: item.sam_file

    cractools:
      brick: CRACTOOLS_EXTRACT
      config:
        SAM_FILE: "{{item.sam_file}}"
        CHIMERA_OUTPUT_PREFIX: "chimeras/cractools/{{item.name}}"
        SPLICE_OUTPUT_PREFIX: "splices/cractools/{{item.name}}"
        MUTATION_OUTPUT_PREFIX: "mutations/cractools/{{item.name}}"
        NB_THREADS: "{{config.nb_threads}}"
      loop: 
        list: samples
        id: name
        export: 
          - value: "{{this.CHIMERA_OUTPUT}}"
            to: tem.chimeras 
          - value: "{{this.SPLICE_OUTPUT}}"
            to: item.splices
          - value: "{{this.MUTATION_OUTPUT}}"
            to: item.mutations

    counts:
      brick: FEATURECOUNTS
      config:
        ANNOTATION_FILE: /data/annotations/human/Homo_sapiens.GRCh38.77.gtf
        INPUT_ALIGNMENTS: "{{samples.extract('sam_file')}}"
        OUTPUT_FILE: "quantification/gene_counts.tsv"

# CREATING A REUSABLE MAKEFILE BRICK - WIP

## PRINCIPLES

If you want to create to reusable bricks, create a new file under brick/
directory with the name of your brick (in capslock by convention) and a
`.Makefile` extension.

In this Makefile template, you need to prefiex the variables names that belong
to the template with `%%`. The only mandatory thing that need to appear is a
rule called `%%`. This is the rule that will be called for this brick in the
generated Makfile.

## THE LATEX EXAMPLE

For exemple, let's create our first re-usable brick that will wrap the
compilation of a latex file. Lets call it `LATEX_PDF.Makefile` :

    %%_LATEX_FILE = undef 
    %%_OUTPUT_FILE = $(basename $(%%_LATEX_FILE)).pdf
    %%_PDFLATEX_OPTIONS:
    %%_PDFLATEX_BINARY: pdflatex

    %%: $(%%_OUTPUT_FILE)

    $(%%_OUTPUT_FILE): $(%%_LATEX_FILE)
      $(%%_PDFLATEX_BINARY) $(%%_PDFLATEX_OPTIONS) $(%%_LATEX_FILE)

This is a very simple wrapper around pdflatex that takes a .tex file and
generate a .pdf from it. There is only one mandatory argument to intanciate
this brick : `LATEX_FILE` wich is set to undef in order to control its
declaration.

As you may know, pdflatex generates a lots of intermediary files during the
compilation. You can implement a clean rule named `%%_clean` that will be added
to the rule `clean` in the finale Makefile. We can add such a rule :

     %%_clean: %%_clean_latex
        rm $(%%_OUTPUT_FILE)

     %%_clean_latex:
       rm -f $(addprefix $(basename $(%%_LATEX_FILE), .aux .bbl .blg .log .out .toc .glsdefs .mtc* .maf .glo .gls .ist .glg .lof))

This should remove the stuff.

One thing we could add to make this brick more usable is a variable that
contains a bibtex file that must be compiled with the file.

    %%_LATEX_FILE = undef
    %%_BIBTEX_FILE =
    %%_OUTPUT_FILE = $(basename $(%%_LATEX_FILE)).pdf
    %%_PDFLATEX_OPTIONS:
    %%_PDFLATEX_BINARY: pdflatex
    %%_BIBTEX_BINARY: bibtex

    %%: $(%%_OUTPUT_FILE)

    $(%%_OUTPUT_FILE): $(%%_LATEX_FILE)
      $(%%_PDFLATEX_BINARY) $(%%_PDFLATEX_OPTIONS) $(%%_LATEX_FILE)
    ifneq ($(strip $(%%_BIBTEX_FILE)),)
      $(%%_BIBTEX_BINARY) $(basename $(%%_LATEX_FILE))
    endif
      $(%%_PDFLATEX_BINARY) $(%%_PDFLATEX_OPTIONS) $(%%_LATEX_FILE)

     %%_clean: %%_clean_latex
        rm $(%%_OUTPUT_FILE)

     %%_clean_latex:
       rm -f $(addprefix $(basename $(%%_LATEX_FILE), .aux .bbl .blg .log .out .toc .glsdefs .mtc* .maf .glo .gls .ist .glg .lof))

We can now use our brand new re-usable brick in a nice pipeline.

    ---
    report:
      brick: LATEX
      config:
        LATEX_FILE: report.tex

## BRICK DOCUMENTATION

You can easily add some documentation to a brick that is parsed and shown to
the user when `./makepipe brick MY_BRICK` is invoked.
Makepipe has no special syntax for documentation, instead it uses the the commentaries
that are placed in the Makefile to build the documentation.

    ## DESCRIPTION
    # The LATEX brick is able to compile a latex file into a pdf using pdflatex command.
    # It is also able to integrate a bibtex file for the bibliography.
    #
    ## VARIABLES
    %%_LATEX_FILE = undef#Latex file to be compiled
    %%_BIBTEX_FILE = #Bibtex file to compile with the LATEX_FILE if needed
    %%_OUTPUT_FILE = $(basename $(%%_LATEX_FILE)).pdf# Output PDF generated by the brick (Cannot be changed)
    %%_PDFLATEX_OPTIONS: #Options to provide to the pdflatex binary
    %%_PDFLATEX_BINARY: pdflatex# pdflatex binary to use
    %%_BIBTEX_BINARY: bibtex# Bibtex binary to use

    %%: $(%%_OUTPUT_FILE)

    $(%%_OUTPUT_FILE): $(%%_LATEX_FILE)
      $(%%_PDFLATEX_BINARY) $(%%_PDFLATEX_OPTIONS) $(%%_LATEX_FILE)
    ifneq ($(strip $(%%_BIBTEX_FILE)),)
      $(%%_BIBTEX_BINARY) $(basename $(%%_LATEX_FILE))
    endif
      $(%%_PDFLATEX_BINARY) $(%%_PDFLATEX_OPTIONS) $(%%_LATEX_FILE)

     %%_clean: %%_clean_latex
        rm $(%%_OUTPUT_FILE)

     %%_clean_latex:
       rm -f $(addprefix $(basename $(%%_LATEX_FILE), .aux .bbl .blg .log .out .toc .glsdefs .mtc* .maf .glo .gls .ist .glg .lof))
