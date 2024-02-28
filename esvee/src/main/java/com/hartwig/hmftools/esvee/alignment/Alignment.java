package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.List;

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.PhaseSet;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class Alignment
{
    private final SvConfig mConfig;

    private final Aligner mAligner;

    public Alignment(final SvConfig config, final Aligner aligner)
    {
        mConfig = config;
        mAligner = aligner;
    }

    public void processPhaseSet(final PhaseSet phaseSet)
    {
        final byte[] fullAssembly = phaseSet.buildFullAssembly();

        if(fullAssembly == null)
            return;

        AssemblyLink assemblyLink = phaseSet.assemblyLinks().get(0);
        List<BwaMemAlignment> alignments = mAligner.alignSequence(fullAssembly);

        for(BwaMemAlignment alignment : alignments)
        {
            SV_LOGGER.debug("assemblyLink({}) alignment: {}", assemblyLink, alignment);
        }

    }

    public static String alignmentStr(final BwaMemAlignment alignment)
    {
        return format("");
    }


}
