from bx.intervals.intersection import Interval,IntervalTree

class CAGEPeak:
    """
    A class to represent and query CAGE (Cap Analysis of Gene Expression) peaks from a BED file.

    Attributes
    ----------
    cage_bed_filename : str
        The filename of the BED file containing CAGE peak data.
    cage_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects containing intervals of peaks.

    Methods
    -------
    __init__(cage_bed_filename):
        Initializes the CAGEPeak object with the given BED filename and reads the BED file to populate the peaks.
    read_bed():
        Reads the BED file and populates the cage_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=10000):
        Queries the CAGE peaks to determine if a given position falls within a peak and calculates the distance to the nearest TSS.
    """
    def __init__(self, cage_bed_filename):
        self._validate_input(cage_bed_filename)
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def _validate_input(self, cage_bed_filename):
        if not cage_bed_filename.endswith('.bed'):
            raise ValueError("CAGE peak file must be in BED format.")
        if not os.path.exists(cage_bed_filename):
            raise FileNotFoundError("CAGE peak file does not exist.")
        
    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split('\t')
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            tss0 = int((start0+end1)/2)
            self.cage_peaks[(chrom,strand)].insert(start0, end1, (tss0, start0, end1))

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within a cage peak>, <nearest distance to TSS>
        If the distance is negative, the query is upstream of the TSS.
        If the query is outside of the peak upstream of it, the distance is NA
        """
        within_peak, dist_peak = 'FALSE', float('inf')
        peaks = self.cage_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for tss0, start0, end1 in peaks:
            # Checks if the TSS is upstream of a peak
            if (strand == '+' and start0 > query and end1 > query) or \
            (strand == '-' and start0 < query and end1 < query):
                continue
            # Checks if the query is within the peak and the distance to the TSS
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            d = (tss0 - query) * (-1 if strand == '-' else 1)
            w = 'TRUE' if within_out else 'FALSE'
            
            if not within_peak=='TRUE':
                within_peak, d = w, (tss0 - query) * (-1 if strand=='-' else +1)
                if within_peak == 'TRUE' or abs(d) < abs(dist_peak):
                    dist_peak = d
                
            else:
                d = (tss0 - query) * (-1 if strand=='-' else +1)
                if abs(d) < abs(dist_peak) and not(w == 'FALSE' and within_peak == 'TRUE'):
                    within_peak, dist_peak = w, d 
        if dist_peak == float('inf'):
            dist_peak = 'NA'
        return within_peak, dist_peak

def read_CAGE_peaks(CAGE_peak):
    if CAGE_peak:
        print("**** Reading CAGE Peak data.", file=sys.stdout)
        return CAGEPeak(CAGE_peak)
    return None

class PolyAPeak:
    """
    A class to represent and query polyA peaks from a BED file.

    Attributes
    ----------
    polya_bed_filename : str
        The filename of the BED file containing polyA peak information.
    polya_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects representing intervals of peaks.

    Methods
    -------
    __init__(polya_bed_filename)
        Initializes the PolyAPeak object with the given BED filename and reads the BED file to populate polyA peaks.
    read_bed()
        Reads the BED file and populates the polya_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=100)
        Queries the polyA peaks to determine if a given position falls within a specified search window of any peak.
    """
    def __init__(self, polya_bed_filename):
        self.polya_bed_filename = polya_bed_filename
        self.polya_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.polya_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            self.polya_peaks[(chrom,strand)].insert(start0, end1, (start0, end1))

    def find(self, chrom, strand, query, search_window=100):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within some distance to polyA>, <distance to closest>
        - if downstream, + if upstream (watch for strand!!!)
        """
        assert strand in ('+', '-')
        within_polyA, dist_polyA = 'FALSE', 'NA'
        hits = self.polya_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for start0, end1 in hits:
            # Checks if the query is within the tail and the distance to the 5'
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            distance = start0 - query if strand == '+' else query - end1

            if within_out:
                within_polyA = 'TRUE'
            if dist_polyA == 'NA' or abs(distance) < abs(dist_polyA):
                dist_polyA = distance

        return within_polyA, dist_polyA

def read_polyA_peaks(polyA_peak):
    if polyA_peak:
        print("**** Reading polyA Peak data.", file=sys.stdout)
        return PolyAPeak(polyA_peak)
    return None