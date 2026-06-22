package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.File;

import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionConfig;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

// Shared read-only inputs handed to every LiftBackWorker: the sidecar resolver and junction index hold no
// mutable state, so all workers share one instance. Each worker builds its own engines + ref-source handle
// from this (RefSequenceSource and the rescue/tail-extend stats are not thread-safe).
public final class LiftBackResources
{
    public final LiftBackResolver Resolver;
    public final AnnotatedJunctionIndex JunctionIndex; // nullable
    public final String RefGenomeFile;
    public final RescueConfig Rescue;
    public final TailExtensionConfig TailExtension;
    public final int TerminalAnchor;
    public final ExcludedRegions ExcludedRegions; // nullable: drop fragments here before lifting

    public LiftBackResources(
            final LiftBackResolver resolver, final AnnotatedJunctionIndex junctionIndex, final String refGenomeFile,
            final RescueConfig rescue, final TailExtensionConfig tailExtension, final int terminalAnchor,
            final ExcludedRegions excludedRegions)
    {
        Resolver = resolver;
        JunctionIndex = junctionIndex;
        RefGenomeFile = refGenomeFile;
        Rescue = rescue;
        TailExtension = tailExtension;
        TerminalAnchor = terminalAnchor;
        ExcludedRegions = excludedRegions;
    }

    // one handle per caller. IndexedFastaSequenceFile is not thread-safe, so workers each open their own.
    public RefSequenceSource openRefSource()
    {
        return openRefSource(RefGenomeFile);
    }

    public static RefSequenceSource openRefSource(final String refGenomeFile)
    {
        if(refGenomeFile == null)
        {
            return null;
        }
        try
        {
            IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(refGenomeFile));
            return (chromosome, posStart, posEnd) ->
            {
                synchronized (fasta)
                {
                    try
                    {
                        ReferenceSequence seq = fasta.getSubsequenceAt(chromosome, posStart, posEnd);
                        return seq != null ? seq.getBases() : null;
                    }
                    catch(Exception e)
                    {
                        return null;
                    }
                }
            };
        }
        catch(Exception e)
        {
            TARS_LOGGER.warn("could not open ref FASTA {} for ref-verify: {}", refGenomeFile, e.toString());
            return null;
        }
    }
}
