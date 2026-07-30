package com.hartwig.hmftools.tars.common;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class TarsConstants
{
    public static final Logger TARS_LOGGER = LogManager.getLogger(TarsConstants.class);

    public static final String APP_NAME = "Tars";

    // suffix appended to a chromosome to name its alt contig, e.g. chr1 -> chr1_tx
    public static final String ALT_CONTIG_SUFFIX = "_tx";

    // bwa-mem scoring, used wherever tars re-scores against the genome
    public static final int MATCH = 1;
    public static final int MISMATCH = -4;
    public static final int GAP_OPEN = -6;
    public static final int GAP_EXTEND = -1;

    // min score to keep a terminal anchor's junction
    public static final int MIN_OVERHANG_SCORE = 5;

    // overhangs longer than this are trusted without scoring
    public static final int MIN_OVERHANG_LENGTH = 12;

    // TODO revisit: OverhangReconciler uses this to suppress genomic scoring for a whole record, which also silences
    // scorable tx-contig candidates; a per-candidate fragment test would remove that use entirely
    public static final int MAX_IMPLIED_INTRON_LENGTH = 1_000_000;
}
