
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - DB Comparision Tool"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.compar.Compar")
}
