
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - DB Comparision Tool"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
}

application {
    mainClass.set("com.hartwig.hmftools.compar.Compar")
}
