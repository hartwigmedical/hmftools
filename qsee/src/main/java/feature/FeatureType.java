package feature;

public enum FeatureType
{
    AMBER_QC                     (VisType.TABLE, SourceTool.AMBER),
    PURPLE_QC                    (VisType.TABLE, SourceTool.PURPLE),
    PURPLE_PURITY                (VisType.TABLE, SourceTool.PURPLE),
    COVERAGE_STATS               (VisType.TABLE, SourceTool.BAM_METRICS),
    READ_STATS                   (VisType.TABLE, SourceTool.BAM_METRICS),

    COVERAGE_DISTRIBUTION        (VisType.PLOT, SourceTool.BAM_METRICS),
    FRAG_LENGTH_DISTRIBUTION     (VisType.PLOT, SourceTool.BAM_METRICS),
    GC_BIAS                      (VisType.PLOT, SourceTool.COBALT),
    DISCORDANT_FRAG_TYPE_COUNTS  (VisType.PLOT, SourceTool.ESVEE),
    DUPLICATE_FREQ               (VisType.PLOT, SourceTool.REDUX),
    MISSED_VARIANT_LIKELIHOOD    (VisType.PLOT, SourceTool.BAM_METRICS),
    BQR_PER_SNV96_CONTEXT        (VisType.PLOT, SourceTool.REDUX),
    BQR_PER_ORIG_QUAL            (VisType.PLOT, SourceTool.REDUX),
    MSI_ERROR_RATES              (VisType.PLOT, SourceTool.REDUX),
    MSI_INDEL_ERROR_BIAS         (VisType.PLOT, SourceTool.REDUX),
    ;

    private final VisType mVisType;
    private final SourceTool mSourceTool;

    FeatureType(VisType visType, SourceTool sourceTool)
    {
        mSourceTool = sourceTool;
        mVisType = visType;
    }

    public SourceTool sourceTool() { return mSourceTool; }
    public VisType visType() { return mVisType; }

    public enum SourceTool
    {
        AMBER,
        BAM_METRICS,
        COBALT,
        ESVEE,
        PURPLE,
        REDUX,
    }

    public enum VisType
    {
        TABLE,
        PLOT,
    }

    @Override
    public String toString()
    {
        return String.format("featureType(%s) sourceTool(%s) visType(%s)",
                this.name(), mSourceTool.name(), mVisType.name());
    }
}
