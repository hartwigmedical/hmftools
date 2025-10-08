package prep;

public class FeatureValue
{
    public String mFeature;
    public String mValue;
    public DataType mDataType;
    public FeatureType mFeatureType;

    public FeatureValue(String feature, String value, DataType dataType, FeatureType featureType)
    {
        mFeature = feature;
        mValue = value;
        mDataType = dataType;
        mFeatureType = featureType;
    }
}
