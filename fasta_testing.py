from classes import AnnotationCreator, SequenceRenderer, FastaSaver

dna_description_fasta = AnnotationCreator.check_format_and_parse('FASTA-1.fasta')
renderer = SequenceRenderer(dna_description_fasta)
renderer.render_forward_strand()
print()
renderer.render_reverse_strand()
print()
renderer.render_double_strand()
print('-' * 500)
dna_description_gb = AnnotationCreator.check_format_and_parse('gene.gb')
renderer = SequenceRenderer(dna_description_gb)
renderer.render_forward_strand()
print()
renderer.render_reverse_strand()
print()
renderer.render_double_strand()

FastaSaver.save_dna_to_fasta('new_fasta.fasta', dna_description_fasta)