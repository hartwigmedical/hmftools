package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.read.ReadContextFactory.createDelContext;
import static com.hartwig.hmftools.sage.read.ReadContextFactory.createInsertContext;
import static com.hartwig.hmftools.sage.read.ReadContextFactory.createSNVContext;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
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

    bcftools annotate -a /Users/jon/hmf/resources/SOMATIC_PON.vcf.gz -c SOMATIC_PON_COUNT ${input}.germline.vcf.gz -O z -o ${input}.pon.vcf.gz
    bcftools index ${input}.pon.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/out_150_hg19.mappability.bed.gz -h /Users/jon/hmf/resources/mappability.hdr -c CHROM,FROM,TO,-,MAPPABILITY ${input}.pon.vcf.gz -O z -o ${input}.map.vcf.gz
    bcftools index ${input}.map.vcf.gz

    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA -c FILTER colo829.sage.map.vcf.gz -O z -o colo829.sage.pre.vcf.gz
    bcftools index colo829.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA -c FILTER colo829.sage.pre.vcf.gz -O z -o colo829.sage.final.vcf.gz


    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA COLO829v003.sage.map.vcf.gz -O z -o COLO829v003.sage.pre.vcf.gz
    bcftools index COLO829v003.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA COLO829v003.sage.pre.vcf.gz -O z -o COLO829v003.sage.final.vcf.gz

    */

    private static final String DISTANCE_FROM_REF_TAG = "NM";
    private final int minQuality;

    private final GenomeRegion bounds;
    private final RefSequence refGenome;
    private final RefContextCandidates candidates;
    private final boolean tumor;
    private final SageConfig config;

    public RefContextConsumer(boolean tumor, final SageConfig config, @NotNull final GenomeRegion bounds,
            @NotNull final RefSequence refGenome, @NotNull final RefContextCandidates candidates) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        this.minQuality = config.minMapQuality();
        this.candidates = candidates;
        this.tumor = tumor;

        this.config = config;
    }

    @Override
    public void accept(final SAMRecord record) {

        if (inBounds(record)) {

            int alignmentStart = record.getAlignmentStart();
            int alignmentEnd = record.getAlignmentEnd();

            if (record.getMappingQuality() >= minQuality) {
                final byte[] refBases = refGenome.alignment(alignmentStart, alignmentEnd);
                final CigarHandler handler = new CigarHandler() {
                    @Override
                    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                            final int refPosition) {
                        processAligned(record, readIndex, refPosition, element.getLength(), refBases);
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
    }

    private void processInsert(@NotNull final CigarElement e, @NotNull final SAMRecord record, int readIndex, int refPosition,
            byte[] refBases) {

        int refIndex = refPosition - record.getAlignmentStart();

        if (refPosition <= bounds.end() && refPosition >= bounds.start()) {
            final String ref = new String(refBases, refIndex, 1);
            final String alt = new String(record.getReadBases(), readIndex, e.getLength() + 1);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null && refContext.readDepth() <= config.maxReadDepth()) {
                if (tumor) {
                    refContext.altRead(ref, alt, createInsertContext(alt, refPosition, readIndex, record, refIndex, refBases));
                } else {
                    refContext.altRead(ref, alt);
                }
            }
        }
    }

    private void processDel(@NotNull final CigarElement e, @NotNull final SAMRecord record, int readIndex, int refPosition,
            byte[] refBases) {

        int refIndex = refPosition - record.getAlignmentStart();

        if (refPosition <= bounds.end() && refPosition >= bounds.start()) {
            final String ref = new String(refBases, refIndex, e.getLength() + 1);
            final String alt = new String(record.getReadBases(), readIndex, 1);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null && refContext.readDepth() <= config.maxReadDepth()) {
                if (tumor) {
                    refContext.altRead(ref, alt, createDelContext(ref, refPosition, readIndex, record, refIndex, refBases));
                } else {
                    refContext.altRead(ref, alt);
                }
            }
        }
    }

    private void processAligned(@NotNull final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength,
            byte[] refBases) {

        int refBasesStartIndex = refPositionStart - record.getAlignmentStart();

        for (int i = 0; i < alignmentLength; i++) {

            int refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refBasesStartIndex + i;

            if (!inBounds(refPosition)) {
                continue;
            }

            final byte refByte = refBases[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];
            final int baseQuality = record.getBaseQualities()[readBaseIndex];

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null && refContext.readDepth() <= config.maxReadDepth()) {

                if (readByte != refByte) {
                    final String alt = String.valueOf((char) readByte);
                    if (tumor) {
                        refContext.altRead(ref, alt, createSNVContext(refPosition, readBaseIndex, record, refBaseIndex, refBases));
                    } else {
                        refContext.altRead(ref, alt);
                    }

                } else {
                    refContext.refRead(record.getMappingQuality(), baseQuality);
                }
            }
        }

    }


    private boolean inBounds(final SAMRecord record) {
        return record.getEnd() >= bounds.start() && record.getStart() <= bounds.end();
    }

    private boolean inBounds(final long position) {
        return position >= bounds.start() && position <= bounds.end();
    }

}
