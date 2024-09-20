import csv
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description="make a tikz plot of the model weights")
parser.add_argument("input", help="the .csv file containing the model weights")
parser.add_argument("--minlinewidth", help="the minimum line width", default=0.3, type=float)
parser.add_argument("--maxlinewidth", help="the maximum line width", default=5.0, type=float)

args = parser.parse_args()



def read_data(filename):
    # Read the data from the CSV file

    with open(filename) as f:
        data = csv.reader(f)
        data = [x for x in data]

    header = data[0]
    data = [row for row in data if row != header]

    index_weight = header.index("weight")
    index_population = header.index("population")
    index_organ = header.index("organ")
    try:
        index_ki = header.index("ki")
    except ValueError as e:
        index_ki = None
    try:
        index_precursor = header.index("precursor")
    except ValueError as e:
        index_precursor = None

    data_dict = {}

    for row in data:
        key = (row[index_population], row[index_organ])
        if key not in data_dict:
            data_dict[key] = {}
        src_key = row[0]
        if index_ki is not None:
            src_key = src_key + row[index_ki]
        if index_precursor is not None:
            src_key = row[index_precursor] + src_key
        data_dict[key][src_key.lower()] = float(row[index_weight])

    return data_dict


# latex preamble


preamble_template = r"""
\documentclass[crop,tikz]{standalone}
\usepackage[english]{babel}
\usepackage{amsmath,amssymb}
\usepackage{bbm}
\usepackage{graphicx}
\usepackage{tikz}
\usepackage{xcolor}
\usepackage{caption}
\usepackage{subcaption}

\usetikzlibrary{arrows, positioning, shapes, fit}

\definecolor{tabblue}{HTML}{1f77b4}
\definecolor{taborange}{HTML}{ff7f0e}
\definecolor{tabgreen}{HTML}{2ca02c}
\definecolor{tabred}{HTML}{d62728}
\definecolor{tabpurple}{HTML}{9467bd}
\definecolor{tabbrown}{HTML}{8c564b}
\definecolor{tabpink}{HTML}{e377c2}
\definecolor{tabgray}{HTML}{7f7f7f}
\definecolor{tabolive}{HTML}{bcbd22}
\definecolor{tabcyan}{HTML}{17becf}


\definecolor{kihi}{HTML}{ffaa00}
\definecolor{kimed}{HTML}{ff0000}
\definecolor{kilo}{HTML}{000000}

\tikzset{
	%Define standard arrow tip
	>=stealth',
	% Define arrow style
	pil/.style={
		-stealth,
		thick,
		shorten <=2pt,
		shorten >=2pt
	},
	map/.style={
		|->,
		thick,
		shorten <=2pt,
		shorten >=2pt
	},
	box/.style = {draw,inner sep=5pt,rounded corners=5pt},
	wbox/.style = {black, above, midway, sloped, yshift=-2pt},
	% legend for mode of entry
	kihinode/.style={shape=rectangle, fill=kihi},
	kimednode/.style={shape=rectangle, fill=kimed},
	kilonode/.style={shape=rectangle, fill=kilo},
	% source node styles
	popnode/.style={circle, draw, black, outer color=white, inner sep=1},
	srcnode/.style={popnode, minimum size=10mm},
	srcnodeN/.style={srcnode, inner color=taborange!50},
	srcnodeCM/.style={srcnode, inner color=tabgreen!50},
	srcnodeEM/.style={srcnode, inner color=tabred!50},
	srcnodeRMneg/.style={srcnode, inner color=tabpurple!50},
	trgtnode/.style={popnode, minimum size=15mm},
	trgtnodeEM/.style={trgtnode, inner color=tabred!50},
	trgtnodeRM/.style={trgtnode, inner color=tabblue!50},
	% edge styles
	edgekihi/.style={color=kihi, bend left=25},
	edgekimed/.style={color=kimed, bend left=0},
	edgekilo/.style={color=kilo, bend left=-25},
}

\newcommand{\ewf}{<<<MAXLINEWIDTH>>>} % edge width factor
\newcommand{\minw}{<<<MINLINEWIDTH>>>} % line width offset

% default model weights

\newcommand{\wRMh}{0.0}
\newcommand{\wNh}{0.0}
\newcommand{\wCMh}{0.0}
\newcommand{\wEMh}{0.0}

\newcommand{\wRMm}{0.0}
\newcommand{\wNm}{0.0}
\newcommand{\wCMm}{0.0}
\newcommand{\wEMm}{0.0}

\newcommand{\wRMl}{0.0}
\newcommand{\wNl}{0.0}
\newcommand{\wCMl}{0.0}
\newcommand{\wEMl}{0.0}


% position ofset of each of the panels

\newcommand{\xloc}{0.0}
\newcommand{\yloc}{0.0}

\begin{document}

\begin{tikzpicture}
% legend for mode of entry

\matrix [below left] at (2.7, -1) {
	\node [kihinode,label=right:{\small Ki67 high}] {}; \\
	\node [kimednode,label=right:{\small Ki67 int}] {}; \\
	\node [kilonode,label=right:{\small Ki67 low}] {}; \\
};

% lewgend for the line widths

\matrix [below left, column sep = 10pt] at (2.7, -7) {
	\node (n11) {}; & \node [label=right:{\small $w = 0.0$}] (n12) {}; \\
	\node (n21) {}; & \node [label=right:{\small $w = 0.5$}] (n22) {}; \\
	\node (n31) {}; & \node [label=right:{\small $w = 1.0$}] (n32) {}; \\
};

\path (n11) edge[line width = \minw+\ewf*0] (n12);
\path (n21) edge[line width = \minw+\ewf*0.5] (n22);
\path (n31) edge[line width = \minw+\ewf*1.0] (n32); 
"""

# insert definition of the model weights

def define_model_weights_EM(dd):
    model_weights_EM = f"""

    \\renewcommand{{\\wNh}}{{{dd["4naivekihi"]:0.2f}}}
    \\renewcommand{{\\wCMh}}{{{dd["4cmkihi"]:0.2f}}}
    \\renewcommand{{\\wEMh}}{{{dd["4emkihi"]:0.2f}}}

    \\renewcommand{{\\wNm}}{{{dd["4naivekimid"]:0.2f}}}
    \\renewcommand{{\\wCMm}}{{{dd["4cmkimid"]:0.2f}}}
    \\renewcommand{{\\wEMm}}{{{dd["4emkimid"]:0.2f}}}

    \\renewcommand{{\\wNl}}{{{dd["4naivekilo"]:0.2f}}}
    \\renewcommand{{\\wCMl}}{{{dd["4cmkilo"]:0.2f}}}
    \\renewcommand{{\\wEMl}}{{{dd["4emkilo"]:0.2f}}}
    """
    return model_weights_EM


def define_model_weights_RM(dd):
    model_weights_RM = f"""

    \\renewcommand{{\\wRMh}}{{{dd["69kihi"]:0.2f}}}
    \\renewcommand{{\\wNh}}{{{dd["4naivekihi"]:0.2f}}}
    \\renewcommand{{\\wCMh}}{{{dd["4cmkihi"]:0.2f}}}
    \\renewcommand{{\\wEMh}}{{{dd["4emkihi"]:0.2f}}}

    \\renewcommand{{\\wRMm}}{{{dd["69kimid"]:0.2f}}}
    \\renewcommand{{\\wNm}}{{{dd["4naivekimid"]:0.2f}}}
    \\renewcommand{{\\wCMm}}{{{dd["4cmkimid"]:0.2f}}}
    \\renewcommand{{\\wEMm}}{{{dd["4emkimid"]:0.2f}}}

    \\renewcommand{{\\wRMl}}{{{dd["69kilo"]:0.2f}}}
    \\renewcommand{{\\wNl}}{{{dd["4naivekilo"]:0.2f}}}
    \\renewcommand{{\\wCMl}}{{{dd["4cmkilo"]:0.2f}}}
    \\renewcommand{{\\wEMl}}{{{dd["4emkilo"]:0.2f}}}
    """
    return model_weights_RM


template_EM = r"""
\node[trgtnodeEM] at (0+\xloc,0+\yloc) (trgt-EM) {<<<LABEL>>>};
\node[srcnodeN] at (-2+\xloc,2+\yloc) (N) {N};
\node[srcnodeCM] at (2+\xloc,2+\yloc) (CM) {CM};
\node[srcnodeEM] at (-2+\xloc,-2+\yloc) (EM) {EM};

\path[pil] (N) edge[edgekihi, line width=\minw+\ewf*\wNh] node[wbox] {\small\wNh} (trgt-EM);
\path[pil] (CM) edge[edgekihi, line width=\minw+\ewf*\wCMh] node[wbox] {\small\wCMh} (trgt-EM);
\path[pil] (EM) edge[edgekihi, line width=\minw+\ewf*\wEMh] node[wbox] {\small\wEMh} (trgt-EM);

\path[pil] (N) edge[edgekimed, line width=\minw+\ewf*\wNm] node[wbox] {\small\wNm} (trgt-EM);
\path[pil] (CM) edge[edgekimed, line width=\minw+\ewf*\wCMm] node[wbox] {\small\wCMm} (trgt-EM);
\path[pil] (EM) edge[edgekimed, line width=\minw+\ewf*\wEMm] node[wbox] {\small\wEMm} (trgt-EM);

\path[pil] (N) edge[edgekilo, line width=\minw+\ewf*\wNl] node[wbox] {\small\wNl} (trgt-EM);
\path[pil] (CM) edge[edgekilo, line width=\minw+\ewf*\wCMl] node[wbox] {\small\wCMl} (trgt-EM);
\path[pil] (EM) edge[edgekilo, line width=\minw+\ewf*\wEMl] node[wbox] {\small\wEMl} (trgt-EM);

\node[box, fit=(N)(CM)(EM)(trgt-EM), label=above:{<<<LABEL>>>}] {};
"""

template_RM = r"""
\node[trgtnodeRM] at (0+\xloc,0+\yloc) (trgt-RM) {<<<LABEL>>>};
\node[srcnodeRMneg] at (2+\xloc,-1.9+\yloc) (CD69neg) {CD69\textsuperscript{-}};
\node[srcnodeN] at (-2+\xloc,2+\yloc) (N) {N};
\node[srcnodeCM] at (2+\xloc,2+\yloc) (CM) {CM};
\node[srcnodeEM] at (-2+\xloc,-2+\yloc) (EM) {EM};

\path[pil] (CD69neg) edge[edgekihi, line width=\minw+\ewf*\wRMh] node[wbox] {\small\wRMh} (trgt-RM);
\path[pil] (N) edge[edgekihi, line width=\minw+\ewf*\wNh] node[wbox] {\small\wNh} (trgt-RM);
\path[pil] (CM) edge[edgekihi, line width=\minw+\ewf*\wCMh] node[wbox] {\small\wCMh} (trgt-RM);
\path[pil] (EM) edge[edgekihi, line width=\minw+\ewf*\wEMh] node[wbox] {\small\wEMh} (trgt-RM);

\path[pil] (CD69neg) edge[edgekimed, line width=\minw+\ewf*\wRMm] node[wbox] {\small\wRMm} (trgt-RM);
\path[pil] (N) edge[edgekimed, line width=\minw+\ewf*\wNm] node[wbox] {\small\wNm} (trgt-RM);
\path[pil] (CM) edge[edgekimed, line width=\minw+\ewf*\wCMm] node[wbox] {\small\wCMm} (trgt-RM);
\path[pil] (EM) edge[edgekimed, line width=\minw+\ewf*\wEMm] node[wbox] {\small\wEMm} (trgt-RM);

\path[pil] (CD69neg) edge[edgekilo, line width=\minw+\ewf*\wRMl] node[wbox] {\small\wRMl} (trgt-RM);
\path[pil] (N) edge[edgekilo, line width=\minw+\ewf*\wNl] node[wbox] {\small\wNl} (trgt-RM);
\path[pil] (CM) edge[edgekilo, line width=\minw+\ewf*\wCMl] node[wbox] {\small\wCMl} (trgt-RM);
\path[pil] (EM) edge[edgekilo, line width=\minw+\ewf*\wEMl] node[wbox] {\small\wEMl} (trgt-RM);


\node[box, fit=(trgt-RM)(CD69neg)(N)(CM)(CM)(EM), label=above:{<<<LABEL>>>}] {};
"""



# skin EM plot

offset_skin_EM = r"""
\renewcommand{\xloc}{0.0}
\renewcommand{\yloc}{0.0}
"""

skin_EM = offset_skin_EM + template_EM.replace("<<<LABEL>>>", "Skin EM")


# skin RM plot

offset_skin_RM = r"""
\renewcommand{\xloc}{6.0}
\renewcommand{\yloc}{0.0}

"""

skin_RM = offset_skin_RM + template_RM.replace("<<<LABEL>>>", "Skin RM")


# LP EM plot 

offset_LP_EM = r"""
\renewcommand{\xloc}{0.0}
\renewcommand{\yloc}{-6.0}

"""

LP_EM = offset_LP_EM + template_EM.replace("<<<LABEL>>>", "LP EM")


# LP RM plot

offset_LP_RM = r"""
\renewcommand{\xloc}{6.0}
\renewcommand{\yloc}{-6.0}

"""

LP_RM = offset_LP_RM + template_RM.replace("<<<LABEL>>>", "LP RM")


# postamble

postamble = r"""
\end{tikzpicture}
\end{document}
"""

def main():
    # import the data
    data_dict = read_data(args.input)

    # define the model weights
    model_weights_skin_RM = define_model_weights_RM(data_dict[("4EM.CD69+", "SK")])
    model_weights_skin_EM = define_model_weights_EM(data_dict[("4EM", "SK")])
    model_weights_LP_RM = define_model_weights_RM(data_dict[("4EM.CD69+", "LP")])
    model_weights_LP_EM = define_model_weights_EM(data_dict[("4EM", "LP")])

    # write the latex code to a file
    texfile = os.path.splitext(os.path.basename(args.input))[0] + ".tex"

    maxlw = str(args.maxlinewidth)
    minlw = str(args.minlinewidth)
    preamble = preamble_template.replace("<<<MAXLINEWIDTH>>>", maxlw).replace("<<<MINLINEWIDTH>>>", minlw)

    with open(texfile, 'w') as f:
        f.write(preamble)
        f.write(model_weights_skin_EM)
        f.write(skin_EM)
        f.write(model_weights_skin_RM)
        f.write(skin_RM)
        f.write(model_weights_LP_EM)
        f.write(LP_EM)
        f.write(model_weights_LP_RM)
        f.write(LP_RM)
        f.write(postamble)

    # compile the latex code
    subprocess.call(['pdflatex', texfile])

if __name__ == "__main__":
    main()