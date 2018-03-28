# Genomics Immerge Trends Analisys (GITA)
  
This project take account for new approaches to analyze genomic data. 


## Installing GITA




### Required

By following those steps, you'll install the application development environment

1. Clone Git repository:
  ```bash
  $ git clone git@github.com:rmf4/GITA.git
  ```
2. The application uses [SAMtools](http://samtools.sourceforge.net/) as main alignment tool.
  ```bash
   $ sudo apt-get install samtools 
  ```

3. Create a [`virtualenv`](https://virtualenv.pypa.io/en/latest/index.html) to host the application:
  You may need `sudo` to install `virtualenv` globally
  ```bash
  # install virtualenv tool manager via pip
  $ [sudo] pip install virtualenv
  # create a new virtualenv folder called venv
  $ virtualenv venv
  # activate your virtualenv!
  $ source venv/bin/activate
  ```

4. Install application dependencies via pip:
  **/!\ Be sure to have your virtualenv activated /!\\**
  This is stipulated by `(venv)` in front of your terminal prompt.

  ```bash
  (venv) $ pip install -r requirements.txt
  ```
5. Usage:

  ```bash
  (venv) $ python trends_for_cnv.py -i1 <control_bam_file> -i2 <case_bam_file> -r <genome_reference_file> -w <window_factor_size> -gff <annotation_file>
  (venv) $ python trends_for_cnv.py -i1 <control_bam_file> -i2 <case_bam_file> -p <ploidy> -w <window_size> -r <genome_reference_file>

  ```
