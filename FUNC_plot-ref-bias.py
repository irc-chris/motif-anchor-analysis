
if __name__ == '__main__':
    print("--------- In Set Plot Info MAIN ---------")

def ishawnia_determine_track(plot_info, track_totals):
    plot_info['track_coords'] = (track_totals[0], track_totals[1])

def ishawnia_ref_alt_expected(strand, phase, plot_info, hic_totals):
    """
    - Assumes that hic_totals is ordered hom1, hom2 contact counts, not dependent of ref or alt
    - strand is the better log, so if REF means score for reference strand (using phase) is better than for alternate
    - phase = 0|1 means reference strand is first so is hom1, 1|0 means reference is second, so is hom2
    - Assumes phase is onlye 0|1 or 0|1
    """

    if phase == "0|1":
        count_r = hic_totals[0]
        count_a = hic_totals[1]
    else:
        count_r = hic_totals[1]
        count_a = hic_totals[0]

    ref_better = count_r >= count_a

    if strand == 'REF' and ref_better:
        plot_info['ra_expect'] = True
    elif strand == 'REF':
        plot_info['ra_expect'] = False
    elif strand == 'ALT' and ref_better:
        plot_info['ra_expect'] = False
    else:
        plot_info['ra_expect'] = True

    plot_info['ra_coords'] = (count_r, count_a)

    if plot_info['ra_expect']:
        return "Homolog with BETTER motif also had MORE contacts (Expected)."
    else:
        return "Homolog with BETTER motif did not also have MORE contacts (Unexpected)."


def ishawnia_determine_expected(plot_info, score_h1, score_h2, hic_h1, hic_h2, one_dim_score, strand):

    # Store whether H1 is the larger for each dimension
    h1_1D_greater = h1_2D_greater = h1_1D_better = motif_expected = False

    score_h1 = float(score_h1)
    score_h2 = float(score_h2)

    h1_1D_greater = score_h1 > score_h2
    h1_2D_greater = hic_h1 > hic_h2
    
    ## Account for if 1D scoring system prefers smaller values
    if one_dim_score == '<':
        print("FLIPPING SCORE DIRECTION")
        h1_1D_better = not h1_1D_greater
    else:
        h1_1D_better = h1_1D_greater

    # If ref performs better on both or worse on both, that's expected
    motif_expected = h1_1D_better == h1_2D_greater
    if motif_expected:
        if h1_1D_better:
            motif_status = "Homolog with BETTER motif has MORE contacts (Expected)"
            plot_info['1D-hmlg'] = 'h1'
            plot_info['2D-hmlg'] = 'h1'
        else:
            motif_status = "Homolog with WORSE motif has FEWER contacts (Expected)"
            plot_info['1D-hmlg'] = 'h2'
            plot_info['2D-hmlg'] = 'h2'
        plot_info['expected'] = True
    else:
        if h1_1D_better:
            motif_status = "Homolog with BETTER motif has FEWER contacts (Unexpected)"
            plot_info['1D-hmlg'] = 'h1'
            plot_info['2D-hmlg'] = 'h2'
        else:
            motif_status = "Homolog with WORSE motif has MORE contacts (Unexpected)"
            plot_info['1D-hmlg'] = 'h2'
            plot_info['2D-hmlg'] = 'h1'
        plot_info['expected'] = False
    
    plot_info['coords'] = (hic_h1, hic_h2)
    plot_info['scores'] = (score_h1, score_h2)

    if strand == "REF" and h1_1D_better:
        plot_info["valid-phase"] = "0|1"
    elif strand == "REF" and not h1_1D_better:
        plot_info["valid-phase"] = "1|0"
    elif strand == "ALT" and h1_1D_better:
        plot_info["valid-phase"] = "1|0"
    else:
        plot_info["valid-phase"] = "0|1"

    return motif_status

def ishawnia_ref_alt_text_colors(strand, phase):
    if strand == "ALT":
        if phase == "0|1":
            text_colors = [0,2,0] # 2 is green; gets better with SNP on PWM
        elif phase == "1|0":
            text_colors = [2,0,0]
        else:
            raise ValueError(f"REF: Invalid phase: {phase}. It must be '1|0' or '0|1'.")
    elif strand == "REF":
        if phase == "0|1":
            text_colors = [0,1,0] # 1 is red; gets worse with SNP on PWM
        elif phase == "1|0":
            text_colors = [1,0,0]
        else:
            raise ValueError(f"ALT: Invalid phase: {phase}. It must be '1|0' or '0|1'.")
    else:
        return ValueError(f"strand {strand} is invalid; strand can only be REF or ALT")
        return text_colors