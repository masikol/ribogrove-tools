
from functools import partial


_organism_name_template_en = '<i>{genus}</i> <i>{species}</i> strain {strain}'
_organism_name_template_ru = '<i>{genus}</i> <i>{species}</i> {strain}'
_organism_name_template_ua = '<i>{genus}</i> <i>{species}</i> {strain}'
_organism_name_template_be = '<i>{genus}</i> <i>{species}</i> {strain}'


def _retrieve_strain_name(strain_name, organism_name_template):
    strain_name = _preprocess_strain_name(strain_name)
    words = strain_name.split(' ')

    is_candidatus = words[0] == 'Candidatus'
    genus_offset = 2 if is_candidatus else 1
    genus_name = _make_genus_name(words, is_candidatus)

    is_sp = 'sp.' in words
    if is_sp:
        sp_name = ' '.join( words[genus_offset+1 :] )
        return '<i>{}</i> sp. {}'.format(genus_name, sp_name)
    else:
        species_name = '{}'.format(words[genus_offset])
    # end if

    species_offset = genus_offset + 1
    strain_designation = ' '.join(
        words[species_offset :]
    )

    return organism_name_template.format(
        genus=genus_name,
        species=species_name,
        strain=strain_designation
    )
# end def


retrieve_strain_name_en = partial(
    _retrieve_strain_name,
    organism_name_template=_organism_name_template_en
)
retrieve_strain_name_ru = partial(
    _retrieve_strain_name,
    organism_name_template=_organism_name_template_ru
)
retrieve_strain_name_ua = partial(
    _retrieve_strain_name,
    organism_name_template=_organism_name_template_ua
)
retrieve_strain_name_be = partial(
    _retrieve_strain_name,
    organism_name_template=_organism_name_template_be
)


def _preprocess_strain_name(strain_name):
    strain_name = _truncate_strain_name(strain_name)
    strain_name = _remove_strain_word(strain_name)
    strain_name = _remove_odd_spaces(strain_name)
    return strain_name
# end def

def _truncate_strain_name(strain_name):
    truncation_marks = [
        ', complete',
        'complete',
        'chromosome',
        'plasmid'
    ]
    for mark in truncation_marks:
        strain_name = strain_name.partition(mark)[0]
    # end for
    return strain_name
# end def

def _remove_strain_word(strain_name):
    return strain_name.replace('strain', ' ')
# end def


def _remove_odd_spaces(strain_name):
    strain_name = strain_name.strip()
    while '  ' in strain_name:
        strain_name = strain_name.replace('  ', ' ')
    # end while
    return strain_name
# end def


def _make_genus_name(strain_name_words, is_candidatus):
    if is_candidatus:
        genus = ' '.join( (strain_name_words[0], strain_name_words[1]) )
        genus = genus.replace('Candidatus', 'Ca.')
    else:
        genus = strain_name_words[0]
    # end if

    return '{}'.format(genus)
# end def
