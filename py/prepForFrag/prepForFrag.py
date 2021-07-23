from filterAnnotation import filterAnnotation

class prepForFrag:
    def __init__(self, l_p_annotated_files, pr_out):
        self.l_p_annotated_filtered_files = l_p_annotated_files
        self.pr_out = pr_out

    def filteredAnnotation(self, p_criteria):

        l_filtered_annotation = []
        for p_annotated_file in self.l_p_annotated_filtered_files:
            print(p_annotated_file)
            c_filteria = filterAnnotation.filterAnnotation(p_annotated_file, p_criteria, self.pr_out)
            c_filteria.loadCriteria()
            p_filtered = c_filteria.filterByCriteriaA()
            l_filtered_annotation.append(p_filtered)
        self.l_filtered_annotations = l_filtered_annotation   
        
    
    def mergeAnnotationForFrag(self):

        return 