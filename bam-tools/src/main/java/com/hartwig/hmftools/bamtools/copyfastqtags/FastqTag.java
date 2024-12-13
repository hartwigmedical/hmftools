package com.hartwig.hmftools.bamtools.copyfastqtags;

import static java.lang.String.format;

import static htsjdk.samtools.util.StringUtil.hexStringToBytes;

public class FastqTag
{
    public final String TagName;
    public final Object Value;

    private FastqTag(final String tagName, final Object value)
    {
        TagName = tagName;
        Value = value;
    }

    public static FastqTag fromFormattedAttribute(final String formattedAttribute)
    {
        String[] components = formattedAttribute.split("[:]");
        if(components.length != 3)
            throw new IllegalArgumentException(format("Formatted attribute is malformed (expected three components when split on ':'): %s", formattedAttribute));

        String name = components[0];
        String type = components[1];
        String strValue = components[2];

        if(type.equals("Z"))
        {
            return new FastqTag(name, strValue);
        }

        if(type.equals("A"))
        {
            Character value = strValue.charAt(0);
            return new FastqTag(name, value);
        }

        if(type.equals("i"))
        {
            Integer value = Integer.parseInt(strValue);
            return new FastqTag(name, value);
        }

        if(type.equals("f"))
        {
            Float value = Float.parseFloat(strValue);
            return new FastqTag(name, value);
        }

        if(type.equals("H"))
        {
            byte[] value = hexStringToBytes(strValue);
            return new FastqTag(name, value);
        }

        throw new IllegalStateException(format("Tag type(%s) is not recognised in attributeString(%s)", type, formattedAttribute));
    }
}
