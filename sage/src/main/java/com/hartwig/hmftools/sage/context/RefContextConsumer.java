package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextFactory;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.CigarHandler;
import com.hartwig.hmftools.sage.sam.CigarTraversal;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RefContextConsumer implements Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(RefContextConsumer.class);

    /*

    bcftools filter -i 'FORMAT/QUAL[1:0]>99 && FORMAT/AD[0:1] < 4' GIABvsSELFv004.sage.vcf.gz -O z -o GIABvsSELFv004.sage.filtered.vcf.gz
    bcftools filter -i 'FORMAT/QUAL[1:0]>99 && FORMAT/AD[0:1] < 4' COLO829v003.sage.vcf.gz -O z -o COLO829v003.sage.filtered.vcf.gz



    bgzip colo829.sage.vcf
    bcftools index colo829.sage.vcf.gz

    input=$1

    bcftools annotate -a /Users/jon/hmf/resources/GERMLINE_PON.vcf.gz -c GERMLINE_PON_COUNT ${input}.vcf.gz -O z -o ${input}.germline.vcf.gz
    bcftools index ${input}.germline.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/SOMATIC_PON.vcf.gz -c SOMATIC_PON_COUNT ${input}.germline.vcf.gz -O z -o ${input}.somatic.vcf.gz
    bcftools index ${input}.somatic.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/SageGermlinePon.hg19.vcf.gz -c PON_COUNT ${input}.somatic.vcf.gz -O z -o ${input}.pon.vcf.gz
    bcftools index ${input}.pon.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/out_150_hg19.mappability.bed.gz -h /Users/jon/hmf/resources/mappability.hdr -c CHROM,FROM,TO,-,MAPPABILITY ${input}.pon.vcf.gz -O z -o ${input}.map.vcf.gz
    bcftools index ${input}.map.vcf.gz


    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA -c FILTER colo829.sage.map.vcf.gz -O z -o colo829.sage.pre.vcf.gz
    bcftools index colo829.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA -c FILTER colo829.sage.pre.vcf.gz -O z -o colo829.sage.final.vcf.gz


    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA COLO829v003.sage.map.vcf.gz -O z -o COLO829v003.sage.pre.vcf.gz
    bcftools index COLO829v003.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA COLO829v003.sage.pre.vcf.gz -O z -o COLO829v003.sage.final.vcf.gz


    bcftools filter -e 'PON_COUNT!= "." && (MIN(PON_COUNT) > 9 || (MIN(PON_COUNT) > 2 && INFO/TIER!="HOTSPOT"))' -s SAGE_PON -m+ COLO829v003.sage.final.vcf.gz -O v -o COLO829v003.sage.filtered.vcf
    bcftools filter -e 'PON_COUNT!= "." && (MIN(PON_COUNT) > 9 || (MIN(PON_COUNT) > 2 && INFO/TIER!="HOTSPOT"))' -s SAGE_PON -m+ GIABvsSELFv004.sage.map.vcf.gz -O v -o GIABvsSELFv004.sage.filtered.vcf

    */

    private final SageConfig config;
    private final GenomeRegion bounds;
    private final RefSequence refGenome;
    private final RefContextFactory candidates;
    private final ReadContextFactory readContextFactory;

    public RefContextConsumer(@NotNull final SageConfig config, @NotNull final GenomeRegion bounds, @NotNull final RefSequence refGenome,
            @NotNull final RefContextFactory candidates) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        this.candidates = candidates;
        this.readContextFactory = new ReadContextFactory(config.readContextFlankSize());

        this.config = config;
    }

    @Override
    public void accept(@NotNull final SAMRecord record) {

        if (inBounds(record) && !reachedDepthLimit(record)) {

            final IndexedBases refBases = refGenome.alignment();

            final CigarHandler handler = new CigarHandler() {
                @Override
                public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                        final int refPosition) {
                    processAlignment(record, readIndex, refPosition, element.getLength(), refBases);
                }

                @Override
                public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                        final int refPosition) {
                    processInsert(element, record, readIndex, refPosition, refBases);
                }

                @Override
                public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                        final int refPosition) {
                    processDel(element, record, readIndex, refPosition, refBases);
                }
            };
            CigarTraversal.traverseCigar(record, handler);
        }
    }

    private void processInsert(@NotNull final CigarElement e, @NotNull final SAMRecord record, int readIndex, int refPosition,
            final IndexedBases refBases) {
        int refIndex = refBases.index(refPosition);

        if (refPosition <= bounds.end() && refPosition >= bounds.start()) {
            final String ref = new String(refBases.bases(), refIndex, 1);
            final String alt = new String(record.getReadBases(), readIndex, e.getLength() + 1);
            boolean findReadContext = findReadContext(readIndex, record);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (!refContext.reachedLimit()) {
                final int baseQuality = baseQuality(readIndex, record, alt.length());
                final ReadContext readContext =
                        findReadContext ? readContextFactory.createInsertContext(alt, refPosition, readIndex, record, refBases) : null;
                refContext.altRead(ref, alt, baseQuality, readContext);
            }
        }
    }

    private void processDel(@NotNull final CigarElement e, @NotNull final SAMRecord record, int readIndex, int refPosition,
            final IndexedBases refBases) {
        int refIndex = refBases.index(refPosition);

        if (refPosition <= bounds.end() && refPosition >= bounds.start()) {
            final String ref = new String(refBases.bases(), refIndex, e.getLength() + 1);
            final String alt = new String(record.getReadBases(), readIndex, 1);
            boolean findReadContext = findReadContext(readIndex, record);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (!refContext.reachedLimit()) {
                final int baseQuality = baseQuality(readIndex, record, 2);
                final ReadContext readContext =
                        findReadContext ? readContextFactory.createDelContext(ref, refPosition, readIndex, record, refBases) : null;
                refContext.altRead(ref, alt, baseQuality, readContext);
            }
        }
    }

    private void processAlignment(@NotNull final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength,
            final IndexedBases refBases) {

        int refIndex = refBases.index(refPositionStart);

        for (int i = 0; i < alignmentLength; i++) {

            int refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refIndex + i;

            if (!inBounds(refPosition)) {
                continue;
            }

            final byte refByte = refBases.bases()[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];
            boolean findReadContext = findReadContext(readBaseIndex, record);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (!refContext.reachedLimit()) {

                int baseQuality = record.getBaseQualities()[readBaseIndex];
                if (readByte != refByte) {
                    final String alt = String.valueOf((char) readByte);
                    final ReadContext readContext =
                            findReadContext ? readContextFactory.createSNVContext(refPosition, readBaseIndex, record, refBases) : null;
                    refContext.altRead(ref, alt, baseQuality, readContext);

                    if (config.mnvEnabled()) {
                        int mnvMaxLength = mnvLength(refPosition,
                                refPositionStart + alignmentLength - 1,
                                readBaseIndex,
                                refBaseIndex,
                                record.getReadBases(),
                                refBases.bases());
                        for (int mnvLength = 2; mnvLength <= mnvMaxLength; mnvLength++) {

                            final String mnvRef = new String(refBases.bases(), refBaseIndex, mnvLength);
                            final String mnvAlt = new String(record.getReadBases(), readBaseIndex, mnvLength);

                            // Only check last base because some subsets may not be valid,
                            // ie CA > TA is not a valid subset of CAC > TAT
                            if (mnvRef.charAt(mnvLength - 1) != mnvAlt.charAt(mnvLength - 1)) {
                                final ReadContext mnvReadContext = findReadContext ? readContextFactory.createMNVContext(refPosition,
                                        readBaseIndex,
                                        mnvLength,
                                        record,
                                        refBases) : null;

                                refContext.altRead(mnvRef, mnvAlt, baseQuality, mnvReadContext);
                            }
                        }
                    }
                } else {
                    refContext.refRead();
                }
            }
        }
    }

    private boolean findReadContext(int readIndex, @NotNull final SAMRecord record) {
        return readIndex >= config.readContextFlankSize() && readIndex < record.getReadLength() - config.readContextFlankSize();
    }

    private int mnvLength(int alignment, int alignmentEnd, int readIndex, int refIndex, byte[] readBases, byte[] refBases) {
        int gap = 0;
        for (int i = 1; alignment + i < alignmentEnd; i++) {
            byte refByte = refBases[refIndex + i];
            byte readByte = readBases[readIndex + i];
            if (refByte == readByte) {
                if (gap == 1) {
                    return i - 1;
                } else {
                    gap++;
                }
            } else {
                gap = 0;
            }
        }

        return alignmentEnd - alignment + 1;
    }

    private int baseQuality(int readIndex, SAMRecord record, int length) {
        int maxIndex = Math.min(readIndex + length, record.getBaseQualities().length) - 1;
        int quality = Integer.MAX_VALUE;
        for (int i = readIndex; i <= maxIndex; i++) {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }

    private boolean inBounds(final SAMRecord record) {
        return record.getEnd() >= bounds.start() && record.getStart() <= bounds.end();
    }

    private boolean inBounds(final long position) {
        return position >= bounds.start() && position <= bounds.end();
    }

    private boolean reachedDepthLimit(@NotNull final SAMRecord record) {
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        RefContext startRefContext = candidates.refContext(bounds.chromosome(), alignmentStart);
        RefContext endRefContext = candidates.refContext(bounds.chromosome(), alignmentEnd);

        return startRefContext.reachedLimit() && endRefContext.reachedLimit();
    }

}
