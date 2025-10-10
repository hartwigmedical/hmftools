package prep;

public class FeatureValue<T>
{
    public String mName;
    public T mValue;
    public FeatureType mType;

    public FeatureValue(String name, T value, FeatureType type)
    {
        mName = name;
        mValue = value;
        mType = type;
    }

    public Class<?> getDataType() {
        return mValue.getClass();
    }
}
