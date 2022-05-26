from abc import ABC, abstractmethod

FEATURE_METADATA_TABULATION = '                     '
HYDROGEN_BOND_SIGN = '|'
FASTA_EXTENSION = '.fasta'
GENBANK_EXTENSION = '.gb'


class AnnotationCreator:
    @staticmethod
    def check_format_and_parse(file_path: str):
        """Check file format and conduct appropriate parsing"""
        if file_path.endswith(FASTA_EXTENSION):
            fasta_file = FastaParser().parse_file_content(file_path)
            return fasta_file
        elif file_path.endswith(GENBANK_EXTENSION):
            gb_file = GenBankParser().parse_file_content(file_path)
            return gb_file
        else:
            print('Unsupported format. Possible formats: .fasta or .gb')


class AbstractParser(ABC):
    @abstractmethod
    def parse_file_content(self, file_path: str):
        pass


class FastaParser(AbstractParser):
    def parse_file_content(self, file_path: str):
        """Parse FASTA file content and create AnnotationFasta object"""
        lines = []
        with open(file_path, 'r') as fasta_file:
            for line in fasta_file:
                line = line.replace('\n', '')
                lines.append(line)
        description_dna = lines[0].split(' ', maxsplit=1)
        dna_identifier = description_dna[0].replace('>', '')
        dna_size = description_dna[1].replace('(', '').replace(')', '')
        forward_strand = ''.join(lines[1:])
        return AnnotationFasta(forward_strand, dna_identifier, dna_size)


class GenBankParser(AbstractParser):
    def parse_file_content(self, file_path: str):
        """Parse GenBank file content and create AnnotationGenBank object"""
        with open(file_path, 'r') as gb_file:
            file_lines = []
            for line in gb_file:  # parse file to get only FEATURES lines
                file_lines.append(line)
            title_list = file_lines[0].split()

        locus = title_list[0]
        sequence_identifier = title_list[1]
        sequence_size = ' '.join(title_list[2:4])
        type_of_sequence = title_list[4]
        topology = title_list[5]
        gen_bank_division = title_list[6]
        modification_date = title_list[7]
        dna_features = self.describe_dna_features(file_lines)
        dna_sequence = self.describe_dna_sequence(file_lines)
        return AnnotationGenBank(dna_sequence, sequence_identifier, sequence_size, locus,
                                 type_of_sequence, topology, gen_bank_division, modification_date, dna_features)

    def describe_dna_features(self, file_lines):
        """Parse feature part of GenBank file and create list of features"""
        list_of_features = []
        start = None
        end = None
        for line in file_lines:
            if 'FEATURES' in line:  # it starts in line which contains FEATURES
                start = file_lines.index(line) + 1
            if 'ORIGIN' in line:  # it ends in line which contains ORIGIN
                end = file_lines.index(line)
        if start is not None and end is not None:
            feature_lines = [e for e in file_lines[start:end]]
            feature_name = ''
            metadata_for_feature = {}
            for line in feature_lines:
                if line.startswith(FEATURE_METADATA_TABULATION):
                    if '=' not in line:  # it means than current line is continue of previous feature
                        previous_value = metadata_for_feature[list(metadata_for_feature)[-1]]
                        current_value = f'{previous_value}{line.strip()}'
                        metadata_for_feature[list(metadata_for_feature)[-1]] = current_value.replace('\n', '').replace("<", "").replace('"', '')
                    else:  # it means than it's a new feature
                        k, v = line.split('=')
                        k = k.replace('/', '').strip()
                        v = v.replace('\n', '').replace("<", "").replace('"', '')
                        metadata_for_feature[k] = v
                else:
                    if feature_name != '':  # it means than it's first iteration, when feature isn't parsed yet
                        feature = Feature(feature_name, metadata_for_feature)
                        list_of_features.append(feature)
                        metadata_for_feature = {}  # 'nulling' metadata for new feature
                    feature_first_line = line.strip().split(' ')
                    feature_name = feature_first_line[0]
                    place = feature_first_line[-1]
                    metadata_for_feature['place'] = place
            feature = Feature(feature_name, metadata_for_feature)
            list_of_features.append(feature)
        return list_of_features

    def describe_dna_sequence(self, file_lines):
        """Parse dna sequence of GenBank file"""

        for line in file_lines:
            if 'ORIGIN' in line:
                start = file_lines.index(line) + 1
        tmp_list_of_dna_sequence = [e for e in file_lines[start:]]
        forward_strand = ''.join(tmp_list_of_dna_sequence)
        forward_strand = ''.join([i.upper() for i in forward_strand if not i.isdigit()]).replace('\n', '').replace(' ', '')
        return forward_strand


class Annotation(ABC):
    @abstractmethod
    def __init__(self, forward_strand_dna: str, sequence_identifier: str, sequence_size: str):
        self.size = sequence_size
        self.identifier = sequence_identifier
        self.forward_strand_dna = forward_strand_dna

    def __str__(self):
        return str(f'forward_strand_dna = {self.forward_strand_dna},size = {self.size}, identifier = {self.identifier}')


class AnnotationFasta(Annotation):
    def __init__(self, forward_strand_dna: str, sequence_identifier: str, sequence_size: str):
        Annotation.__init__(self, forward_strand_dna, sequence_identifier, sequence_size)


class AnnotationGenBank(Annotation):
    def __init__(self, forward_strand_dna, sequence_identifier, sequence_size, locus, type_of_sequence,
                 topology, gen_bank_division, modification_date, features):
        Annotation.__init__(self, forward_strand_dna, sequence_identifier, sequence_size)
        self.locus = locus
        self.type_of_sequence = type_of_sequence
        self.topology = topology
        self.gen_bank_division = gen_bank_division
        self.modification_date = modification_date
        self.features = features


class Feature:
    def __init__(self, name: str, metadata: dict):
        """Creates feature object with appropriate metadata"""
        self.name = name
        self.metadata = metadata


class SequenceRenderer:
    def __init__(self, dna_description: AnnotationFasta or AnnotationGenBank):
        self.dna_description = dna_description

    def create_forward_strand(self):
        return self.dna_description.forward_strand_dna

    def render_forward_strand(self):
        """Prints to console forward strand DNA"""
        print(self.create_forward_strand())

    def create_reverse_strand(self):
        """Constructs reverse strand DNA in accordance with the principle of complementarity """
        reverse_strand_nucleotides = []
        for nucleotide in self.dna_description.forward_strand_dna:
            if nucleotide == 'A':
                reverse_strand_nucleotides.append('T')
            elif nucleotide == 'T':
                reverse_strand_nucleotides.append('A')
            elif nucleotide == 'G':
                reverse_strand_nucleotides.append('C')
            elif nucleotide == 'C':
                reverse_strand_nucleotides.append('G')
        reverse_strand_dna = ''.join(reverse_strand_nucleotides)
        return reverse_strand_dna

    def render_reverse_strand(self):
        """Prints to console reverse strand DNA"""
        print(self.create_reverse_strand())

    def create_double_strand(self):
        """Binds DNA complementary chains(forward and reverse)"""
        reverse_and_forward_strands = []
        dna_reverse_strand = self.create_reverse_strand()
        reverse_and_forward_strands.append(dna_reverse_strand)
        reverse_and_forward_strands.append(self.dna_description.forward_strand_dna)
        sign = '\n' + (HYDROGEN_BOND_SIGN * len(self.dna_description.forward_strand_dna)) + '\n'
        double_strand_dna = sign.join(reverse_and_forward_strands)
        return double_strand_dna

    def render_double_strand(self):
        """Prints to console double strand DNA"""
        print(self.create_double_strand())


class FastaSaver:

    @staticmethod
    def save_dna_to_fasta(file_path: str, dna: AnnotationFasta):
        """Saves DNA and sequence description to provided file"""
        with open(file_path, 'a') as file:
            file.write(f'>{dna.identifier} ({dna.size})\n{dna.forward_strand_dna}')



