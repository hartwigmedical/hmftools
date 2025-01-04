package com.hartwig.hmftools.geneutils.mapping;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.geneutils.mapping.SequenceTester.FLD_SEQUENCE;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class SequenceInfo
{
    public final int Id;
    public final String Sequence;
    public final ChrBaseRegion Region;

    public SequenceInfo(final int id, final String sequence, final ChrBaseRegion region)
    {
        Id = id;
        Sequence = sequence;
        Region = region;
    }

    public String toString()
    {
        return format("id(%d) region(%ss) seq(%d:%s)", Id, Region, Sequence.length(), Sequence);
    }

    public static void addSequenceHeader(final StringJoiner sj)
    {
        sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END);
        sj.add(FLD_SEQUENCE);
    }

    public void addSequenceInfo(final StringJoiner sj)
    {
        if(Region != null)
        {
            sj.add(Region.Chromosome);
            sj.add(String.valueOf(Region.start()));
            sj.add(String.valueOf(Region.end()));
        }
        else
        {
            sj.add(null);
            sj.add(null);
            sj.add(null);
        }

        sj.add(Sequence);
    }
}
