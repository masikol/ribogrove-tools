
import re
from functools import partial


_organism_name_template_en = '<i>{genus}</i> <i>{species}</i> strain {strain}'
_organism_name_template_ru = '<i>{genus}</i> <i>{species}</i> {strain}'
_organism_name_template_ua = '<i>{genus}</i> <i>{species}</i> {strain}'
_organism_name_template_be = '<i>{genus}</i> <i>{species}</i> {strain}'


def _retrieve_strain_name(strain_name, asm_acc, organism_name_template):
    strain_name = _preprocess_strain_name(strain_name)
    words = strain_name.split(' ')

    is_candidatus = words[0] == 'Candidatus'
    genus_offset = 2 if is_candidatus else 1
    genus_name = _make_genus_name(words, is_candidatus)

    is_sp = 'sp.' in words
    if is_sp:
        sp_name = ' '.join( words[genus_offset+1 :] )
        raw_strain_name = '<i>{}</i> sp. {}'.format(genus_name, sp_name)
        return _amend_strain_name(
            raw_strain_name,
            asm_acc,
            organism_name_template
        )
    else:
        species_name = '{}'.format(words[genus_offset])
    # end if

    species_offset = genus_offset + 1
    strain_designation = ' '.join(
        words[species_offset :]
    )

    raw_strain_name = organism_name_template.format(
        genus=genus_name,
        species=species_name,
        strain=strain_designation
    )

    return _amend_strain_name(raw_strain_name, asm_acc, organism_name_template)
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


def italicize_candidatus(tax_name):
    # For e.g. 'Candidatus Saccharibacteria'
    regex_match = re.search(r'Ca\. ([A-Z][a-z]+)', tax_name)
    if not regex_match is None:
        return tax_name.replace(
            regex_match.group(0),
            '<i>Ca.</i> {}'.format(regex_match.group(1))
        )
    else:
        return '<i>{}</i>'.format(tax_name)
    # end if
# end def




def _amend_strain_name(strain_name, asm_acc, organism_name_template):
    # Fix known troubles with strain names

    strain_name = _deitalicize_candidatus_gs(strain_name)

    if asm_acc == 'GCF_020541245.1':
        return _deitalicize_candidatus_gs(
            organism_name_template.format(
                genus='Ca. Annandia',
                species='pinicola',
                strain='Ad13-065'
            )
        )
    elif asm_acc == 'GCF_902387285.1':
        return organism_name_template.format(
            genus='Methanomethylophilus',
            species='alvi',
            strain='MGYG-HGUT-02456'
        )
    elif asm_acc == 'GCF_902384015.1':
        return organism_name_template.format(
            genus='Methanosphaera',
            species='stadtmanae',
            strain='MGYG-HGUT-02164'
        )
    elif asm_acc == 'GCF_963668415.1':
        return 'Uncultured <i>Methanospirillum</i> sp.'
    elif asm_acc == 'GCF_963668475.1':
        return 'Uncultured <i>Methanospirillum</i> sp.'
    elif asm_acc == 'GCF_963674485.1':
        return 'Uncultured <i>Methanolobus</i> sp.'
    elif asm_acc == 'GCF_009217585.1':
        return '{} <sup>**</sup>'.format(strain_name)
    elif 'ATCC 33223 ATCC 33223; 39E' in strain_name:
        return strain_name.replace(
            'ATCC 33223 ATCC 33223; 39E',
            'ATCC 33223'
        )
    elif 'SY1 = JCM 12091 SY1' in strain_name:
        return strain_name.replace(
            'SY1 = JCM 12091 SY1',
            'SY1'
        )
    elif 'str. ' in strain_name:
        return strain_name.replace(
            'str. ',
            ''
        )
    elif 'pv. oryzae ' in strain_name:
        return strain_name.replace(
            'pv. oryzae ',
            ''
        )
    # end if

    subsp_regex_match = re.search(r'subsp. [a-z]+ ', strain_name)
    if not subsp_regex_match is None:
        subsp_needle = subsp_regex_match.group(0)
        return strain_name.replace(subsp_needle, '')
    # end if

    return strain_name
# end def


def _deitalicize_candidatus_gs(strain_name):
    # gs means 'genus and species'
    # For e.g. 'Ca. Mancarchaeum acidiphilum strain Mia14'
    candidatus_re_match = re.search(
        r'^\<i\>Ca\. ([A-Z][a-z]+)\</i\> \<i\>([a-z]+)\</i\>',
        strain_name
    )
    if not candidatus_re_match is None:
        genus_name   = candidatus_re_match.group(1)
        species_name = candidatus_re_match.group(2)
        strain_name = strain_name.replace(
            candidatus_re_match.group(0),
            '<i>Ca.</i> {} {}'.format(genus_name, species_name)
        )
    # end if

    # For e.g.
    #          'Ca. Micrarchaeum sp. A_DKE'
    candidatus_re_match = re.search(
        r'^\<i\>Ca\. ([A-Z][a-z]+)\</i\> sp.',
        strain_name
    )
    if not candidatus_re_match is None:
        genus_name   = candidatus_re_match.group(1)
        strain_name = strain_name.replace(
            candidatus_re_match.group(0),
            '<i>Ca.</i> {} sp.'.format(genus_name)
        )
    # end if

    return strain_name
# end if


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
