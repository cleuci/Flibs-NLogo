breed   [flibs flib]     ;; FLiBs (finite living blobs) are the agents of the model: they are structured as
                         ;; finite automata (1, 2). They have "a finite number of states; an input signal
                         ;; causes it to change automatically from one state to another.
                         ;; The kind of automaton used in a flib also generates signals. Incoming and outgoing
                         ;; signals are represented within the automaton by symbols. When a signal is received,
                         ;; the automaton changes its state and emits a second signal" (2). A flib is considered
                         ;; perfect predictor i.e. very well adapted when its outgoing signals (previsions)
                         ;; is equal to the next incoming environmental signal. Signals are always binary digits.


flibs-own [chromosome    ;; every flib owns a chromosome that is a string coding its state transition table;
           state         ;; the current state of the flib
           fitness]      ;; a measure of its forestalling ability


globals [counter         ;; a counter useful in different procedures
         new-input       ;; the incoming environmental signal
         average         ;; mean of flibs population performances
         best            ;; the best flibs performance value
         worst           ;; the worst flibs performance value
         donor           ;; a performing flib sharing part of its genes
         recipient       ;; a flib acquiring genes from a donor
         a-split         ;; a first cut in the recipient's chromosome
         b-split         ;; a second cut in the recipient's chromosome
         wild-type            ;; a natural state chromosome
         optimal]        ;; the flib's chromosome showing a perfect prevision ability



;;  ----------   SETUP PROCEDURES   -----------------------------------------------------------------------------
;;  -------------------------------------------------------------------------------------------------------------


to setup                 ;; initializing the model
  clear-all
  reset-ticks
  set counter 0
  create-flibs num-flibs [
    set shape "face happy"
    set size 1.8
    set chromosome ""
    set label fitness
    fd 2 + random 7]
  ask flibs [chromosome-genesis]
end


to chromosome-genesis   ;; chromosomes are strings randomly built through an iterative procedure: even string's
                        ;; positions are 0 or 1: they represent a possible outgoing signal; odd string's
                        ;; positions represent possible flib's states determined by NUM-STATES slider
  set chromosome word chromosome (word (random 2) (random num-states))
  set counter counter + 1
  ifelse counter < num-states * 2 [chromosome-genesis]
    [set counter 0]
end



;;  ----------   RUNTIME PROCEDURES   ---------------------------------------------------------------------------
;;  -------------------------------------------------------------------------------------------------------------

to search
  ;; checking chromosome dimension
  ask flibs [if length chromosome != num-states * 4
    [show "error: invalid chromosome lenght" stop]]

  ;; testing flibs' prevision ability along 100 cycles and results are analyzed
  set counter 0
  ask flibs [set fitness 0 set state 0]
  ask flibs [chrom-test]
  performances-analysis

  ;; a new generation begins
  tick
  output-print ""
  output-print word "tick: " ticks

  ;; if flibs population harbours a perfect predictor, the goal of the genetic algorithm is reached
  if best = 100 [final-message stop]

  ;; if flibs population doesn't harbour any perfect predictor, genetic algorithm activates two operators:
  ;; genetic shuffling (as a consequence of mating processes) and mutagenesis

  ;; every generation one conjugation event occurs by default, this frequency can be
  ;; lowered by slider "MATE-RATE"
  ifelse random-float 1 < mate-rate [conjugation]
    [output-print "no mating event"]

  ;; every generation, one mutagenesis process occurs  by default, this frequency can be
  ;; lowered by slider MUTATION-RATE
  ifelse random-float 1 < mutation-rate
    [set wild-type [who] of one-of flibs ask flib wild-type [mutate]]
    [output-print "no mutation event"]
end

to final-message
  output-print "        EUREKA!!!"
  set optimal [who] of one-of flibs with [fitness = 100]
  ifelse best > 9  [output-type word "id: " optimal] [output-type word "id: 0" optimal]
  output-type word " optimal predictor:         " [chromosome] of flib optimal
  output-print word "     fitness: " [fitness] of flib optimal
  output-print ""
  output-print "--- no genetic variation (see to previous generation) ---"
  output-print ""
end


;; Operator 1: FITNESS EVALUATION
;; --------------------------------------------------------------------------------------------------------------

to chrom-test
  ;; an environmental input symbol is read from the given sequence
  set new-input read-from-string item (counter mod (length environmental-cycle)) environmental-cycle
  ;; only 0 or 1 symbols are accepted
  if new-input != 0 and new-input != 1 [
    show "error: invalid environmental value"
    stop]
  ;; each flibs prevision ability is tested and fitness is updated
  if read-from-string  item (4 * state + 2 * new-input) chromosome =
        read-from-string item ((counter + 1) mod (length environmental-cycle)) environmental-cycle
            [set fitness (fitness + 1)]
  ;; new flibs state is computed
  set state read-from-string item (4 * state + 2 * new-input + 1) chromosome
  set counter counter + 1
  ;; the cycle is repeated 100 times
  if counter < 100 [chrom-test]
  set counter 0
end

to performances-analysis
  ask flibs [
  set average mean [fitness] of flibs
  set best max [fitness] of flibs
  set worst min [fitness] of flibs
  set label fitness
  ;; flibs performances are visually represented through emoticons language
  ifelse fitness < 55
    [set shape "face sad"] [ifelse fitness < 85
    [set shape "face neutral"]
    [set shape "face happy"]]
  ]
end


;; Operator 2: GENETIC SHUFFLING
;; --------------------------------------------------------------------------------------------------------------

to conjugation
  clear-links
  ;; the conjugation process requires two flibs' chromosomes: the donor and the recipient
  ;; just the second one undergoes to hybridization
  select-flibs
  ;; a link highlighted the two mating flibs
  ask flib donor [create-link-to flib recipient]
  ask flib recipient [genetic-shuffling]
end

to select-flibs
  set donor [who] of one-of flibs with [fitness = best]
  set recipient [who] of one-of flibs
  if donor = recipient [select-flibs]     ;; self conjugation is forbidden
  ;; generating report about the selection process
  ifelse donor > 9  [output-type word "id: " donor] [output-type word "id: 0" donor]
  output-type word " donor chromosome:      " [chromosome] of flib donor
  output-print word "     fitness: " [fitness] of flib donor
  ifelse recipient > 9  [output-type word "id: " recipient]  [output-type word "id: 0" recipient]
  output-type word " recipient chromosome:  " [chromosome] of flib recipient
  output-print word "     fitness: " [fitness] of flib recipient
end

to genetic-shuffling
  ;; a gene's sequence included between a-split and b-split restriction sites is randomly chosen
  set a-split random (num-states * 4)
  set b-split random (num-states * 4)
  if a-split = b-split [genetic-shuffling]
  set counter 0
  hybridization
  ;; generating a report about genetic shuffling
  ifelse recipient > 9  [output-type word "id: " recipient]  [output-type word "id: 0" recipient]
  output-print word " hybridized chromosome: " chromosome
  output-type "- flanking positions of replaced fragment:  [ "
  output-type word a-split " - "
  output-print word b-split " [ "
end

to hybridization
 ;; the genes' sequence included between a-split and b-split restriction sites on a random
 ;; flib is replaced by the corresponding sequence on one of the most performing flibs;
 ;; chromosomes are treated as circle-shaped
 if a-split < b-split [
    set chromosome replace-item (a-split + counter)
    chromosome (item (a-split + counter) [chromosome] of flib donor)
    set counter (counter + 1) if counter < b-split - a-split
    [hybridization]]
 if a-split > b-split [
    set chromosome replace-item ((a-split + counter) mod (num-states * 4))
    chromosome (item ((a-split + counter) mod (num-states * 4)) [chromosome] of flib donor)
    set counter (counter + 1) if counter < (num-states * 4) - (a-split - b-split)
    [hybridization]]
  if length chromosome != num-states * 4 [show "Invalid chromosome length"]
  set counter 0
end


;; Operator 3: MUTAGENESIS
;; --------------------------------------------------------------------------------------------------------------

to mutate
  ;; managing mutation and generating a report
  ifelse wild-type > 9  [output-type word "id: " wild-type]  [output-type word "id: 0" wild-type]
  output-type word " wild-type chromosome:  " [chromosome] of flib wild-type
  output-print word "     fitness: " [fitness] of flib wild-type
  ;; mutations occur randomly at a given frequency on just one locus
  let dice random length chromosome
  let muton read-from-string item dice chromosome
  ifelse dice mod 2 = 0
    [set chromosome replace-item dice chromosome word ((muton + 1) mod 2) ""]
    [set chromosome replace-item dice chromosome word ((muton + 1) mod num-states) ""]
  ifelse wild-type > 9  [output-type word "id: " wild-type] [output-type word "id: 0" wild-type]
  output-print word " mutated chromosome:    " [chromosome] of flib wild-type
end
@#$#@#$#@
GRAPHICS-WINDOW
239
11
520
293
-1
-1
13.0
1
10
1
1
1
0
0
0
1
-10
10
-10
10
0
0
1
ticks
30.0

BUTTON
17
236
111
269
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
127
217
219
250
NIL
search
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
19
143
220
203
environmental-cycle
10010110
1
0
String

TEXTBOX
21
96
214
142
Insert or modify the environmental cycle below. Warning: the only symbols allowed are 0 and 1 digits.
11
104.0
0

SLIDER
17
53
218
86
num-states
num-states
1
10
4.0
1
1
NIL
HORIZONTAL

SLIDER
17
10
218
43
num-flibs
num-flibs
1
100
15.0
1
1
NIL
HORIZONTAL

PLOT
18
306
829
588
performances variations
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"best" 1.0 0 -955883 true "" "if ticks > 0 [plot best]"
"average" 1.0 0 -13791810 true "" "if ticks > 0 [plot average]"
"worst" 1.0 0 -10263788 true "" "if ticks > 0 [plot worst]"

PLOT
540
10
828
291
performances distribution
NIL
NIL
0.0
110.0
0.0
15.0
true
false
"" ""
PENS
"default" 1.0 1 -13840069 true "" "histogram [fitness] of flibs"

BUTTON
127
258
219
291
NIL
search
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
16
674
217
707
mutation-rate
mutation-rate
0
1
0.3
0.01
1
NIL
HORIZONTAL

OUTPUT
240
602
830
707
11

TEXTBOX
48
600
261
634
     Current genetic variations:
13
104.0
1

TEXTBOX
19
605
234
635
_____________________________
12
104.0
1

SLIDER
16
634
217
667
mate-rate
mate-rate
0
1
1.0
0.01
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

Flibs'NLogo reproduces in NetLogo programming language and its integrated modeling environment, a genetic algorithm whose purpose is evolving a perfect predictor from a pool of digital creatures constituted by finite automata (1) or flibs (finite living blobs) that are the agents of the model. The project is based on the rules described by  Alexander K. Dewdney in "Exploring the field of genetic algorithms in a primordial computer sea full of flibs" from the vintage Scientific American column "Computer Recreations" (2, 3) and recollected in the book "The Armchair Universe" (4)
As the same author summarized: "Flibs [...] attempt to predict changes in their environment. In the primordial computer soup, during each generation, the best predictor crosses chromosomes with a randomly selected flib. Increasingly accurate predictors evolve until a perfect one emerges. A flib [...] has a finite number of states, and for each signal it receives (a 0 or a 1) it sends a signal and enters a new state. The signal sent by a flib during each cycle of operation is its prediction of the next signal to be received from the environment" (3)

## HOW IT WORKS

External signals follow cyclically a pattern of binary digits that can be inserted by the user. Even if they could be represented more properly for a human mind by a table or a graph, the repertoires of responses to environmental signals and the state-transitions of flibs are codified linearly in the chromosomal string as well as the genetic information is codified in DNA strand.
The topology of flib's chromosome should be considered circular, a gene is a symbol in the chromosome and an allele is a gene at a particular locus. Flibs live in a kind of primordial soup where time flows discretely. A time unit (tick) represents a generation; each generation is made up of 100 cycles. At every cycle, flibs receive an incoming signal and in dependence of its state compute an outcoming binary signal, entering into a new state, according to the transition table codified by its chromosome. The performances of a flib will be mesured on the basis of the accordance between the outcoming signal and the next incoming environmental signal: this is the fitness operator.
In addition to the fitness operator, a genetic algorithm involves at least two other types of operators: crossover and mutation. Mutation occurs at a given frequency when an allele is changed at random (think of a gene being struck by a cosmic ray); crossover, in this case, is substituted by a more generic genetic shuffling: it occurs through a kind of unidirectional genes exchange when a chromosomic fragment included between two loci randomly chosen on one of the most performing flibs, replaces the corresponding alleles of another flib.

## HOW TO USE IT

The SETUP button initializes the model creating a number of flibs determined by the NUM-FLIBS slider; each flib will be provided with a chromosome representing a casual transition table built on the basis of states number determined by the NUM-STATES slider. 
The SEARCH buttons start the genetic algorithm evaluating first the flibs' fitness and identifying the flibs with relative best and worst performances. If no flib reaches the absolute perfect forestalling performance, the algorithm executes the genetic shuffling between a flib's chromosome with current best fitness and another one randomly chosen and finally it operates the mutation of one random flib's gene. After that, the cycle starts again and a new generation begins.
The model window shows flibs as emoticons labelled with their respective fitness value, a link highlights the mating flibs. There are two plots: the first shows a histogram representing the flibs' fitness distribution; the second one is a three-line graph showing generation after generation, the variation in the best, worst and average flibs' performances. In the bottom, an output window reports details about the last genetic variations. 


## THINGS TO NOTICE

In spite of the model described by Dewdney, Flibs'NLogo exhibits a genetic shuffling different from a true crossing-over: indeed just the receiving chromosome (not always less performing) varies its genetic kit; on the contrary the best flib donates part of its chromosome to the second one, so that no flib is removed from the soup, but one of them changes its structure by hybridization. The mechanism was already used in the model "Minimal Genetic Algorithm".

 
## THINGS TO TRY

Explore the relationship between the length of environmental signals pattern and the flibs' number of states to get the best predictor flib within a given number of generation (or before your patience is exhausted). Obviously, there is a limit to the highest theoretical fitness (the adaptability that in this case means forestalling ability) in a given environmental signals pattern: this theme is discussed by Dewdney (2)
The success of genetic algorithms often depends greatly on some details as the size of the population and the probabilities of genetic shuffling and mutation (5). Test the influence of these parameters acting on the sliders NUM-FLIBS, MATE-RATE, and MUTATION-RATE. Is adaptation easier in large populations where a great number of candidate answers to environmental challenges are available, or in a small population where genetic drift can offer an unexpected chance to adaptation? In flibs' primordial sea, are more influential mutations or genetic shuffling in evolution; in other words: are flibs "hopeful monsters" or gradualist species?

 

## EXTENDING THE MODEL

The genetic variability of flibs population has an influence over the probability to obtain a perfect predictor. An operator exploring this feature could enlighten the relationship between biodiversity and adaptability. 
Dewdney notices: "the concept of a flib is so flexible that input and output can represent a great variety of specific biological phenomena. For example, an input signal could represent a chemical or temperature gradient. The corresponding output symbol could be a command to an effector that controls cilia, or a spore forming mechanism" (2). The model could be extended from a conceptual point of view: some scholars think of biological adaptation as a cognitive process or vice versa, in this case, individual and collective flibs' activity could be regarded as a mental process, as well (6).
 


## NETLOGO FEATURES

Flibs' chromosomes have the double feature as strings and as ordered numeric digits collection. To switch them each other, two NetLogo commands are very useful: “read-from-string” interprets the given chromosome and reports the resulting numerical value; on the contrary “word” command allows to revert a number into a string (because it requires two inputs, the second one could be double quotation marks: "") A third command is useful in mutation and genetic shuffling operators: “replace-item” allows an easy single gene variation as well as a substitution with a chromosome fragment from a donor turtle into the corresponding loci of a recipient turtle.
Besides, "mod" function helps to treat strings (e.g. chromosomes and environmental signals pattern) as circular structures.

## RELATED MODELS

Leuci C. (2018). NetLogo Minimal Genetic Algorithm model http://modelingcommons.org/browse/one_model/5742


## CREDITS AND REFERENCES

1. Brian Hayes "On the finite-state machine, a minimal model of mousetraps, ribosomes and the human soul" Scientific American, Vol. 249, No. 6 (December 1983), pages 19-28. Available at: https://tinyurl.com/ycukjtor
2. A. K. Dewdney “Exploring the field of genetic algorithms in a primordial computer see full of flibs” Scientific American Vol. 253, No. 5 (November 1985), pages 21-33. Partially available at: https://tinyurl.com/y84bnear
3. A. K. Dewdney "The king (a chess program) is dead, long live the king (a chess machine)" Scientific American Vol. 254, No. 2 (February 1986), page 21. Available at: https://tinyurl.com/ycphyqqm
4. A. K. Dewdney “The Armchair Universe: An Exploration of Computer Worlds" 1988 W. H. Freeman & Co. New York
5. Melanie Mitchell "An introduction to genetic algorithms" 1999 MIT Press Ltd, Cambridge, Massachusetts. Available at: http://www.boente.eti.br/fuzzy/ebook-fuzzy-mitchell.pdf
6. Karl Popper "All life is problem solving" 1999 Routledge, New York. Available at: https://tinyurl.com/yc5s2y3e


## COPYRIGHT AND LICENSE

Copyright 2018 Cosimo Leuci.

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
