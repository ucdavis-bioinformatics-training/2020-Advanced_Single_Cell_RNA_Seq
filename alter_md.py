import argparse

def main(infile):
    f = open(infile)
    out = open(infile.split('.')[0] + "_fixed.md", 'w')
    r_found = False
    tick_reset = False
    for line in f:
        # Pass over the r ```r {CONTENT} ```
        if line == "```r\n" and r_found is False:
            r_found = True
            out.write(line)
        elif line == "```\n" and r_found is True:
            r_found = False
            out.write(line)
        # Reset the normal ``` {CONTENT} ```
        elif line == "```\n" and tick_reset is False:
            tick_reset = True
            out.write("<div class='r_output'>")
        elif line == "```\n" and tick_reset is True:
            tick_reset = False
            out.write("</div>")
        else:
            # if we are in an output area get rid of the extra ##
            if tick_reset is True:
                line = line.replace("##", '')
            out.write(line)



parser = argparse.ArgumentParser(description='Alter rmd.py: script to fix rmd files for better rendering.'
                                             'Create an HTML using knitr from the Rmd file and keep_md = True.'
                                             'Then run this script on the file. Output is same name with prefix _fixed',
                                 epilog='For questions or comments, please contact Matt Settles <settles@ucdavis.edu> '
                                        'or Keith Mitchell <kgmitchell@ucdavis.edu\n', add_help=True)
parser.add_argument('-i', '--input', help="Input .md rendered from Rmd to be fixed for better presentation.")

# parser.add_argument('-o', '--output', help="Output .md fixed for better presentation in template.",
#                     type=str)

options = parser.parse_args()
main(options.input)