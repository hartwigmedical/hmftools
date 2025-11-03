package feature;

public enum FeatureType
{
    GENERAL                      (VisType.TABLE),
    MUTATIONAL_BURDEN            (VisType.TABLE),
    COVERAGE_STATS               (VisType.TABLE),
    READ_STATS                   (VisType.TABLE),

    COVERAGE_DISTRIBUTION        (VisType.PLOT),
    FRAG_LENGTH_DISTRIBUTION     (VisType.PLOT),
    GC_BIAS                      (VisType.PLOT),
    DISCORDANT_FRAG_TYPE_COUNTS  (VisType.PLOT),
    DUPLICATE_FREQ               (VisType.PLOT),
    MISSED_VARIANT_LIKELIHOOD    (VisType.PLOT),
    BQR_PER_SNV96_CONTEXT        (VisType.PLOT),
    BQR_PER_ORIG_QUAL            (VisType.PLOT),
    MS_INDEL_ERROR_RATES         (VisType.PLOT),
    MS_INDEL_ERROR_BIAS          (VisType.PLOT),
    ;

    private final VisType mVisType;

    FeatureType(VisType visType)
    {
        mVisType = visType;
    }

    public VisType visType() { return mVisType; }

    public enum VisType
    {
        TABLE,
        PLOT,
    }

    @Override
    public String toString()
    {
        return String.format("featureType(%s) visType(%s)",
                this.name(), mVisType.name());
    }
}
