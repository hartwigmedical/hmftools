package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.CodingUtils.determineCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determineNonCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determinePreOrPostCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determineUpstreamContext;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

public class CodingContext
{
    public TranscriptRegionType RegionType; // favours more impactful type
    public TranscriptCodingType CodingType;
    public int ExonRank;

    public int CodingBase; // indexed at 1, coding base from start of transcript's coding region
    public int[] CodingPositionRange;
    public boolean SpansSpiceJunction;
    public int DeletedCodingBases;

    public int NonCodingBaseDistance;
    public int UpstreamPhase;
    public int BasesToLastExonJunction;

    public CodingContext()
    {
        RegionType = TranscriptRegionType.UNKNOWN;
        CodingType = TranscriptCodingType.UNKNOWN;
        ExonRank = 0;
        SpansSpiceJunction = false;
        CodingPositionRange = new int[] {0, 0};
        DeletedCodingBases = 0;
        NonCodingBaseDistance = 0;
        UpstreamPhase = PHASE_NONE;
        BasesToLastExonJunction = -1;
    }

    public boolean isCoding() { return CodingType == CODING && RegionType == EXONIC; }

    public String hgvsStr() { return "tbc"; }

    public static String csvHeader()
    {
        return "HgvsCoding,RegionType,CodingType,ExonRank,CodingBase,UpstreamPhase,CodingPosRange,SpansSplice,CodingDeleted,NonCodingBaseDist";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(hgvsStr());
        sj.add(RegionType.toString());
        sj.add(CodingType.toString());
        sj.add(String.valueOf(ExonRank));
        sj.add(String.valueOf(CodingBase));
        sj.add(String.valueOf(UpstreamPhase));
        sj.add(CodingPositionRange[SE_START] == CodingPositionRange[SE_END] ? String.valueOf(CodingPositionRange[SE_START])
                : String.format("%d-%d", CodingPositionRange[SE_START], CodingPositionRange[SE_END]));
        sj.add(String.valueOf(SpansSpiceJunction));
        sj.add(String.valueOf(DeletedCodingBases));
        sj.add(String.valueOf(NonCodingBaseDistance));

        return sj.toString();
    }
}
