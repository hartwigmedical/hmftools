package com.hartwig.hmftools.fastqstats;

public class TrackerKey implements Comparable<TrackerKey>{
    private TrackerType type;
    private String name;

    public TrackerKey(TrackerType type, String name){
        this.type = type;
        this.name = name;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        TrackerKey that = (TrackerKey) o;

        return type == that.type && name.equals(that.name);
    }

    @Override
    public int hashCode() {
        return 17 * type.ordinal() + name.hashCode();
    }

    @Override
    public int compareTo(final TrackerKey o) {
        if(type == o.type) {
            return name.compareTo(o.name);
        }
        else {
            return type.compareTo(o.type);
        }
    }

    @Override
    public String toString(){
        return type + " " + name;
    }

}
