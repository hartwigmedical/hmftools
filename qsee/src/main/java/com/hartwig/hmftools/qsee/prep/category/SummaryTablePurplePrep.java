package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.MB_PER_GENOME;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.CONTAMINATION;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.DELETED_GENES;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.LOH_PERCENT;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.PLOIDY;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.PURITY;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.TINC;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.TMB_MS_INDELS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.TMB_SMALL_VARIANTS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.TMB_STRUCTURAL_VARIANTS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.UNSUPPORTED_CN_SEGMENTS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.putFeature;

import java.io.IOException;
import java.util.EnumMap;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;
import com.hartwig.hmftools.qsee.status.PurpleQCStatusConverter;
import com.hartwig.hmftools.qsee.status.QcStatus;
import com.hartwig.hmftools.qsee.status.ThresholdRegistry;

import org.jetbrains.annotations.NotNull;

public class SummaryTablePurplePrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.PURPLE;

    public SummaryTablePurplePrep(QseePrepConfig config) { mConfig = config; }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private PurityContext loadPurplePurity(String sampleId) throws IOException
    {
        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, sampleId);
        String qcFile = PurpleQCFile.generateFilename(baseDir, sampleId);

        return PurityContextFile.readWithQC(qcFile, purityFile);
    }

    @VisibleForTesting
    static EnumMap<PurpleQCStatus, QcStatus> qcStatusFrom(PurityContext purityContext, ThresholdRegistry thresholds)
    {
        EnumMap<PurpleQCStatus, QcStatus> qcStatuses = new EnumMap<>(PurpleQCStatus.class);

        PurpleQCStatusConverter converter = new PurpleQCStatusConverter(thresholds);

        for(PurpleQCStatus purpleQcStatus : PurpleQCStatus.values())
        {
            boolean purpleQcStatusExists = purityContext.qc().status().contains(purpleQcStatus);

            QcStatus qcStatus = purpleQcStatusExists
                    ? converter.toQseeQcStatus(purpleQcStatus)
                    : QcStatus.createEmpty();

            qcStatuses.put(purpleQcStatus, qcStatus);
        }

        return qcStatuses;
    }

    private static QcStatus getTincQcStatus(EnumMap<PurpleQCStatus, QcStatus> qcStatuses)
    {
        QcStatus tincStatus = qcStatuses.get(PurpleQCStatus.FAIL_NO_TUMOR);

        if(tincStatus == null)
            tincStatus = qcStatuses.get(PurpleQCStatus.WARN_TINC);

        if(tincStatus == null)
            tincStatus = QcStatus.createEmpty();

        return tincStatus;
    }

    @VisibleForTesting
    static List<Feature> createFeatures(PurityContext purityContext, ThresholdRegistry thresholds)
    {
        PurpleQC purpleQC = purityContext.qc();
        FittedPurity purpleFit = purityContext.bestFit();

        EnumMap<PurpleQCStatus, QcStatus> qcStatuses = qcStatusFrom(purityContext, thresholds);

        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);

        putFeature(featuresMap, PLOIDY, purpleFit.ploidy(), QcStatus.createEmpty());
        putFeature(featuresMap, PURITY, purpleFit.purity(), qcStatuses.get(PurpleQCStatus.WARN_LOW_PURITY));
        putFeature(featuresMap, LOH_PERCENT, purpleQC.lohPercent(), QcStatus.createEmpty());
        putFeature(featuresMap, UNSUPPORTED_CN_SEGMENTS, purpleQC.unsupportedCopyNumberSegments(), qcStatuses.get(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE));
        putFeature(featuresMap, DELETED_GENES, purpleQC.deletedGenes(), qcStatuses.get(PurpleQCStatus.WARN_DELETED_GENES));

        putFeature(featuresMap, TINC, purpleQC.tincLevel(), getTincQcStatus(qcStatuses));
        putFeature(featuresMap, CONTAMINATION, purpleQC.contamination(), qcStatuses.get(PurpleQCStatus.FAIL_CONTAMINATION));

        putFeature(featuresMap, TMB_SMALL_VARIANTS, purityContext.tumorMutationalBurdenPerMb(), QcStatus.createEmpty());
        putFeature(featuresMap, TMB_MS_INDELS, purityContext.microsatelliteIndelsPerMb(), QcStatus.createEmpty());
        putFeature(featuresMap, TMB_STRUCTURAL_VARIANTS, purityContext.svTumorMutationalBurden() / MB_PER_GENOME, QcStatus.createEmpty());

        return featuresMap.values().stream().toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        if(sampleType != SampleType.TUMOR)
            return List.of();

        PurityContext purityContext = loadPurplePurity(sampleId);
        List<Feature> features = createFeatures(purityContext, mConfig.QcThresholds);
        return features;
    }
}
