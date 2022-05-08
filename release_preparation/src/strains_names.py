
from functools import partial


_organism_name_template_en = '<i>{genus}</i> <i>{species}</i> strain {strain}'
_organism_name_template_ru = '<i>{genus}</i> <i>{species}</i> {strain}'
_organism_name_template_ua = '<i>{genus}</i> <i>{species}</i> {strain}'


def _retrieve_strain_name(title, organism_name_template):
    title = _preprocess_title(title)
    words = title.split(' ')

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
    strain_name = ' '.join( words[species_offset :] )

    return organism_name_template.format(
        genus=genus_name,
        species=species_name,
        strain=strain_name
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


def _preprocess_title(title):
    title = _truncate_title(title)
    title = _remove_strain_word(title)
    title = _remove_odd_spaces(title)
    return title
# end def

def _truncate_title(title):
    truncation_marks = [
        ', complete',
        'complete',
        'chromosome',
        'plasmid'
    ]
    for mark in truncation_marks:
        title = title.partition(mark)[0]
    # end for
    return title
# end def

def _remove_strain_word(title):
    return title.replace('strain', ' ')
# end def


def _remove_odd_spaces(title):
    title = title.strip()
    while '  ' in title:
        title = title.replace('  ', ' ')
    # end while
    return title
# end def


def _make_genus_name(title_words, is_candidatus):
    if is_candidatus:
        genus = ' '.join( (title_words[0], title_words[1]) )
        genus = genus.replace('Candidatus', 'Ca.')
    else:
        genus = title_words[0]
    # end if

    return '<i>{}</i>'.format(genus)
# end def
