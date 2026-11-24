
import polars as pl


def remove_invalid_species(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col('Species').map_elements(
            is_validlike_species_name,
            return_dtype=pl.Boolean
        ).alias('species_valid')
    ).filter(pl.col('species_valid') == True)
# end def


def is_validlike_species_name(input_arg):
    if isinstance(input_arg, str):
        return _seems_like_valid(input_arg)
    elif  isinstance(input_arg, pl.Series):
       return pl.Series([
            _seems_like_valid(x) for x in input_arg
       ])
    else:
        raise TypeError('Incompatible input_arg type: "{}"'.format(type(input_arg)))
   # end if
# end def


def _seems_like_valid(species_str):
    return not species_str.endswith('sp.') \
        and not ' sp. ' in species_str
# end def



# Funciton for top-10 selection

def get_top_rows(domain_df: pl.DataFrame,
                 grouped_df: pl.DataFrame,
                 sorted_traits: list,
                 trait_col_name: str = 'len',
                 top_num: int = 10) -> list[dict]:
    # Collect rows until we have at least top_num genomes, including ties
    top_rows = []
    top_i = 0
    trait_i = 0

    while top_i < top_num:
        target_trait = sorted_traits[trait_i]
        # Get all genomes with this trait
        genomes_at_trait_df = grouped_df.filter(
            pl.col(trait_col_name) == target_trait
        )

        for j in range(genomes_at_trait_df.height):
            row = genomes_at_trait_df.row(j, named=True)
            curr_asm_acc = row['asm_acc']
            curr_trait = row[trait_col_name]

            # Get all rows matching asm_acc and curr_trait
            curr_genome_df = domain_df.filter(
                (pl.col('asm_acc') == curr_asm_acc) & 
                (pl.col(trait_col_name) == curr_trait)
            )

            # Take first row for scalar values, collect seqID list
            first_row = curr_genome_df.row(0, named=True)
            seq_ids = curr_genome_df['seqID'].to_list()

            top_rows.append({
                'asm_acc': first_row['asm_acc'],
                trait_col_name: first_row[trait_col_name],
                'seqID': seq_ids,
                'strain_name': first_row['strain_name'],
                'Domain': first_row['Domain'],
            })
            top_i += 1
        # end for
        trait_i += 1
    # end for

    return top_rows
# end def
