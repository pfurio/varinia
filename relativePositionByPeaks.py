# -*- coding: utf-8 -*-

# Author: Pedro Furió Tarí
# Data: 04/02/2013
#
# Este script cogerá el gtf para saber la posición de cada uno de los genes y los picos detectados por MACS, Pugh ...
# Deolverá 5 ficheros:
#      1. Fichero que contenga la asignación de un pico a un gen
#      2. Vector indicando cuantos picos han sido detectados en cualquier posicion upstream del gen que haya sido
#         asignado relativo a su TSS
#      3. Vector indicando cuantos picos han sido detectados en cualquier posicion del cuerpo del gen relativo al TSS
#      4. Vector indicando cuantos picos han sido detectados en cualquier posicion downstream del gen que haya sido
#         asignado relativo a su TTS
#      5. Vector indicando cuantos picos han sido detectados en cualquier posicion del cuerpo del gen relativo al TTS

# Se asume que el fichero gff tiene este formato:
#chr1    SGD     gene    335     649     .       +       .       ID=YAL069W;Name=YAL069W;Ontology_term=GO:0003674,GO:0005575,GO:0008150;Note=Dubious%20open%20reading%20frame%20unlikely%20to%20encode%20a%20protein%2C%20based%20on%20available%20experimental%20and%20comparative%20sequence%20data;dbxref=SGD:S000002143;orf_classification=Dubious
#chr1    SGD     gene    538     792     .       +       .       ID=YAL068W-A;Name=YAL068W-A;Ontology_term=GO:0003674,GO:0005575,GO:0008150;Note=Dubious%20open%20reading%20frame%20unlikely%20to%20encode%20a%20protein%3B%20identified%20by%20gene-trapping%2C%20microarray-based%20expression%20analysis%2C%20and%20genome-wide%20homology%20searching;dbxref=SGD:S000028594;orf_classification=Dubious

# Los ficheros de picos que se aceptan actualmente pueden ser los generados por Pugh y los que genera MACS2 (narrowPeak)

import getopt, sys, os.path


class Mygenes:
    def __init__(self, gene, start, end, strand):
        self.gene = gene
        self.start = start
        self.end = end
        self.strand = strand

    def getGene(self):
        return self.gene

    def getStart(self):
        return self.start

    def getEnd(self):
        return self.end

    def getStrand(self):
        return self.strand


class MyPeaks(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def getPos(self):
        return [self.start, self.end]


class Pugh(MyPeaks):
    def __init__(self, cw, counts, start, end):
        super(Pugh, self).__init__(start, end)
        self.cw = cw
        self.counts = counts

    def getCW(self):
        return self.cw

    def getCounts(self):
        return self.counts

    def getExtraInfo(self):
        return [self.cw, self.counts]


class MACS(MyPeaks):
    def __init__(self, fold, score, start, end):
        super(MACS, self).__init__(start, end)
        self.score = score
        self.fold = fold

    def getFoldChange(self):
        return self.fold

    def getScore(self):
        return self.score

    def getExtraInfo(self):
        return [self.score, self.fold]


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:p:o:t:b:s", ["help", "gff=", "peaks=", "output=", "type=", "bps=", "score"])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    gff  = None
    peak = None
    out  = None
    typeP = "Pugh"
    bps  = 5000
    score = False

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-g", "--gff"):
            if os.path.isfile(a):
                gff = a
        elif o in ("-p", "--peaks"):
            if os.path.isfile(a):
                peak = a
        elif o in ("-o", "--output"):
            out = a
        elif o in ("-t", "--type"):
            typeP = a
        elif o in ("-b", "--bps"):
            bps = int(a)
        elif o in ("-s", "--score"):
            score = True
        else:
            assert False, "Unhandled option"

    type_opts = ["Pugh","MACS"]

    if gff is not None and out is not None and peak is not None and typeP in type_opts:
        run(gff, peak, out, typeP, bps, score)
    else:
        usage()


def usage():
    print "\nUsage: python relativePositionByPeaks [options] <mandatory>"
    print "Options:"
    print "\t-h, --help:\n\t\t show this help message and exit"
    print "\t-t, --type:\n\t\t Type of peaks given - Pugh or MACS - (default: Pugh)"
    print "\t-b, --bps:\n\t\t Number of base pairs to take into account in order to perform peak-gene pairing"
    print "\t-s, --score:\n\t\t Put this option if you want to count not only the peaks but the counts that fall" \
          " within the peaks"
    print "Mandatory:"
    print "\t-g, --gff:\n\t\t GFF file containing the annotations of the genes"
    print "\t-p, --peaks:\n\t\t Peak file (narrowPeak from MACS or gff from Pugh)"
    print "\t-o, --output:\n\t\t directory to put the output files"
    print "\n04/02/2013. Pedro Furió Tarí.\n"


def parsePughPeaks(gff):

    peaks = {}
    inputPeak = open(gff, 'r')

    contador = 0
    last_chrom = None
    last_start = None
    last_end   = None
    for linea in inputPeak:
        l = linea.split("\t")
        chrom  = l[0]
        start  = int(l[3])
        end    = int(l[4])
        counts = float(l[5])
        cwdistance = int(l[8].split("=")[1])

        if chrom not in peaks:
            peaks[chrom] = []

        if chrom != last_chrom or start != last_start or end != last_end:
            contador += 1
            peaks[chrom].append(Pugh( cwdistance, counts, start, end ))

        last_chrom = chrom
        last_start = start
        last_end   = end

    print str(contador)+ " peaks !!"

    inputPeak.close()

    return peaks


def parseNarrowPeaks(narrowPeak):

    # Only get peaks with score above 100
    peaks = {}
    inputPeak = open(narrowPeak, 'r')

    for linea in inputPeak:
        l = linea.split("\t")
        chrom  = l[0]
        start  = int(l[1])
        end    = int(l[2])
        signal   = float(l[6])
        score  = int(l[4])

        #if score >= 100:
        if chrom not in peaks:
            peaks[chrom] = []

        peaks[chrom].append(MACS(signal, score, start, end ))

    inputPeak.close()

    return peaks


def parseGTF(gff):
    inputGTF = open(gff, 'r')
    genes = {}

    for line in inputGTF:
        linea_split = line.split("\t")
        chrom  = linea_split[0]
        start  = int(linea_split[3])
        end    = int(linea_split[4])
        strand = linea_split[6]

        popurri = linea_split[8]
        gene_id = popurri.split('ID=')[1].split(';')[0]

        # We save the information read
        if chrom not in genes:
            genes[chrom] = []
        genes[chrom].append(Mygenes(gene_id, start, end, strand))

    inputGTF.close()
    return genes


def run(gff, peak, out, typeP, bps, score):

    print "Parsing peak file..."
    # 1. First, we store the peaks
    peaks = None
    if typeP == "Pugh":
        peaks = parsePughPeaks(peak)
    elif typeP == "MACS":
        peaks = parseNarrowPeaks(peak)


    print "Parsing gene file..."
    # 2. Second, we save all the genes with their positions
    genes = parseGTF(gff)

    outSum = open(out + "/PeakGene.txt", 'w')
    # ExtraInfo could be the score and the fold-change for MACS peaks
    # or cw distance and the counts for Pugh peaks

    outSum.write("#chrom\tPeakStart\tPeakEnd\tExtraInfo1\tExtraInfo2\tGeneID\tGeneStart\tGeneEnd\tStrand\tDistance\tPosition\n")
    interg_tss = [0] * bps
    body_tss   = [0] * bps
    body_tts   = [0] * bps
    interg_tts = [0] * bps

    inicio = None
    final  = None

    print "Establishing peak-gene pairs..."
    # 3. Look for peak - gene relations
    for chrom in peaks:

        old_index = 0

        for pico in peaks[chrom]:
            [pstart, pend] = pico.getPos()
            [extra1, extra2] = pico.getExtraInfo()
            extraString = str(extra1) + "\t" + str(extra2)

            # flag indica si el pico cae o no dentro del cuerpo de un gen
            flag = False

            # En lastGene tengo la información del gen más cercano upstream o downstream
            # lastGene = ["gene","strand","start","end"]
            lastGene = ["","",0,0]
            lastGeneFlag = False

            for i in range(old_index, len(genes[chrom])):
                gen = genes[chrom][i]
                strand = gen.getStrand()
                start = gen.getStart()
                end = gen.getEnd()
                genId = gen.getGene()

                #### Possible states of a peak: Being before, inside or after the gene
                #### Each state may also have different substates

                ### Peak before the gene ###
                if pend < start:

                    # Check if the gene before the peak can be also taken into account
                    if lastGeneFlag is True:
                        tipo = "DOWNSTREAM" if lastGene[1] == "+" else "UPSTREAM"

                        if tipo == "UPSTREAM" and pstart - lastGene[3] <=  400 or \
                                                tipo == "DOWNSTREAM" and pstart - lastGene[3] <= 300:

                            if flag is True: # Es válido ... Comprobar lo del flag especial. Está ya en gene_body?
                                tipo += "_GB"

                            outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                         lastGene[0] + "\t" + str(lastGene[2]) + "\t" + str(lastGene[3]) + "\t"
                                         + lastGene[1] + "\t" + str(pstart - lastGene[3]) + "\t" + tipo + "\n")

                            inicio = pstart - lastGene[3]
                            final  = pend - lastGene[3]

                            if lastGene[1] == "-":
                                if score is True:
                                    for w in range(min(inicio,bps), min(final,bps)):
                                        interg_tss[w] += extra2
                                else:
                                    for w in range(min(inicio,bps), min(final,bps)):
                                        interg_tss[w] += 1
                            else:
                                if score is True:
                                    for w in range(min(inicio,bps), min(final,bps)):
                                        interg_tts[w] += extra2
                                else:
                                    for w in range(min(inicio,bps), min(final,bps)):
                                        interg_tts[w] += 1


                    # Check the current gene
                    tipo = "UPSTREAM" if strand == "+" else "DOWNSTREAM"

                    if tipo == "UPSTREAM" and start - pend <=  400 or \
                                            tipo == "DOWNSTREAM" and start - pend <= 300:

                        if flag is True: # Es válido ... Comprobar lo del flag especial. Está ya en gene_body?
                            tipo += "_GB"

                        outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                     genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                     + strand + "\t" + str(start - pend) + "\t" + tipo + "\n")

                        inicio = start - pend
                        final  = start - pstart

                        if strand == "+":
                            if score is True:
                                for w in range(min(inicio,bps), min(final,bps)):
                                    interg_tss[w] += extra2
                            else:
                                for w in range(min(inicio,bps), min(final,bps)):
                                    interg_tss[w] += 1
                        else:
                            if score is True:
                                for w in range(min(inicio,bps), min(final,bps)):
                                    interg_tts[w] += extra2
                            else:
                                for w in range(min(inicio,bps), min(final,bps)):
                                    interg_tts[w] += 1

                    # End loop. We have everything we want from the peak
                    break

                ### Peak in the body ###

                elif pstart >= start and pend <= end:
                    # Peak completely inside the gene

                    flag = True

                    distance_tss = 0
                    distance_tts = 0

                    if strand == "+":

                        # TSS
                        distance_tss = pstart - start

                        inicio = pstart - start
                        final  = pend - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += 1

                        # TTS
                        distance_tts = end - pend

                        inicio = end - pend
                        final  = end - pstart

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += 1

                    else:

                        # TSS
                        distance_tss = end - pend

                        inicio = end - pend
                        final  = end - pstart

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += 1

                        # TTS
                        distance_tts = pstart - start

                        inicio = pstart - start
                        final  = pend - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += 1

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t" +
                                 strand + "\t" + str(distance_tss) + "\tGENE_BODY_TSS\n")

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t" +
                                 strand + "\t" + str(distance_tts) + "\tGENE_BODY_TTS\n")

                elif start >= pstart and end <= pend:
                    # Gene completely inside the peak
                    print genId + " gene completely inside the peak..."

                    flag = True

                    distance_tss = 0
                    distance_tts = 0

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                 + strand + "\t" + str(distance_tss) + "\tGENE_BODY_TSS\n")

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                 + strand + "\t" + str(distance_tts) + "\tGENE_BODY_TTS\n")

                elif start <= pstart <= end:
                    # Peak partially inside the gene

                    flag = True

                    distance_tss = 0
                    distance_tts = 0

                    if strand == "+":
                        distance_tss = pstart - start

                        # Body TSS
                        inicio = pstart - start
                        final  = end - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += 1


                        distance_tts = 0

                        # Body TTS
                        inicio = 0
                        final  = end - pstart

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += 1

                        # Downstream
                        inicio = 0
                        final  = pend - end

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tts[w] += 1

                        outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                     genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                     + strand + "\t0\tDOWNSTREAM\n")

                    else:

                        # Body TSS
                        distance_tss = 0

                        inicio = 0
                        final  = end - pstart

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += 1


                        # Body TTS
                        distance_tts = pstart - start

                        inicio = pstart - start
                        final  = end - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += 1

                        # Upstream
                        inicio = 0
                        final  = pend - end

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tss[w] += 1

                        outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                     genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                     + strand + "\t0\tUPSTREAM\n")

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                 + strand + "\t" + str(distance_tss) + "\tGENE_BODY_TSS\n")

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                 + strand + "\t" + str(distance_tts) + "\tGENE_BODY_TTS\n")

                elif start <= pend <= end:
                    # Peak partially inside the gene
                    #    <----------------->
                    # XXXXXXX

                    flag = True

                    distance_tss = 0
                    distance_tts = 0

                    if strand == "+":
                        distance_tss = 0

                        # Body TSS
                        inicio = 0
                        final  = pend - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += 1


                        distance_tts = end - pend

                        # Body TTS
                        inicio = end - pend
                        final  = end - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += 1

                        # Upstream
                        inicio = 0
                        final  = start - pstart

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tss[w] += 1

                        outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                     genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                     + strand + "\t0\tUPSTREAM\n")

                    else:

                        # Body TSS
                        distance_tss = end - pend

                        inicio = end - pend
                        final  = end - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tss[w] += 1


                        # Body TTS
                        distance_tts = 0

                        inicio = 0
                        final  = pend - start

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                body_tts[w] += 1

                        # Downstream
                        inicio = 0
                        final  = start - pstart

                        if score is True:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tts[w] += extra2
                        else:
                            for w in range(min(inicio,bps), min(final,bps)):
                                interg_tts[w] += 1

                        outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                     genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                     + strand + "\t0\tDOWNSTREAM\n")

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                 + strand + "\t" + str(distance_tss) + "\tGENE_BODY_TSS\n")

                    outSum.write(chrom + "\t" + str(pstart) + "\t" + str(pend) + "\t" + extraString + "\t" +
                                 genId + "\t" + str(start) + "\t" + str(end) + "\t"
                                 + strand + "\t" + str(distance_tts) + "\tGENE_BODY_TTS\n")


                ### Peak after the gene ###
                elif end < pstart:
                    if pstart - bps > start:
                        old_index = i

                    lastGeneFlag = True
                    lastGene = [genId, strand, start, end]

    outSum.close()

    print "Creating plot files..."

    upstr  = open(out + "/upstream.txt", 'w')
    for i in interg_tss:
        upstr.write(str(i) + "\n")
    upstr.close()

    down = open(out + "/downstream.txt", 'w')
    for i in interg_tts:
        down.write(str(i) + "\n")
    down.close()

    bodys   = open(out + "/geneBody_TSS.txt", 'w')
    for i in body_tss:
        bodys.write(str(i) + "\n")
    bodys.close()

    bodys   = open(out + "/geneBody_TTS.txt", 'w')
    for i in body_tts:
        bodys.write(str(i) + "\n")
    bodys.close()



if __name__ == "__main__":
    main()




