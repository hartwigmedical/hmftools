package com.hartwig.hmftools.tars.common;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class TarsConstants
{
    public static final Logger TARS_LOGGER = LogManager.getLogger(TarsConstants.class);

    public static final String APP_NAME = "Tars";

    // alt contig name suffix, e.g. chr1 -> chr1_tx
    public static final String ALT_CONTIG_SUFFIX = "_tx";

    // a deletion this size or smaller either side of a splice N is folded into it (10M100N5D5M -> 10M105N5M) so the
    // junction stays on the annotated boundary; anything larger is kept as a real deletion
    public static final int MAX_MERGED_DELETION_BP = 5;

    // bwa-mem scoring, used wherever tars re-scores against the genome
    public static final int MATCH = 1;
    public static final int MISMATCH = -4;
    public static final int GAP_OPEN = -6;
    public static final int GAP_EXTEND = -1;

    // min score to keep a terminal anchor's junction
    public static final int MIN_OVERHANG_SCORE = 5;

    // overhangs longer than this are trusted without scoring
    public static final int MIN_OVERHANG_LENGTH = 12;

    // MAPQ of a placement taken as confident
    public static final int CONFIDENT_MAPQ = 60;

    // bwa's own default -T floor: a primary still under it after lift-back is unmapped, a supplementary only bwa's
    // lowered -T 19 kept is dropped
    public static final int PRIMARY_AS_UNMAP_THRESHOLD = 30;
    public static final int SUPP_AS_DROP_THRESHOLD = 30;

    // Prefer a tied locus near its mate.
    public static final int MATE_PROXIMITY_MAX_DISTANCE = 1_000_000;

    // Valid intron bounds for supplementary merges.
    public static final int MIN_IMPLIED_INTRON_LENGTH = 21;

    public static final int MAX_IMPLIED_INTRON_LENGTH = 1_000_000;

    // Most supplementaries folded into one primary.
    public static final int MAX_SUPP_MERGES = 2;

    // Maximum read-span overlap accepted during a supplementary merge.
    public static final int MAX_SUPP_READ_OVERLAP = 5;

    // ref-verify only: how far the primary's own boundary may be pulled back into its soft clip while
    // looking for an annotated junction to splice at
    public static final int MAX_REF_VERIFY_BOUNDARY_SHIFT = 8;
}
