"""
Constants for column names and configuration
"""

# Column names
class ColumnName:
    ID = 'ID'
    RANK = 'Rank'
    READS = 'Reads'
    RPU = 'RPU'
    SEQUENCES = 'sequences'
    LENGTH='length'
    DISTANCE = 'Distance'

    # cluster
    CLUSTER = 'Cluster'
    RANK_IN_CLUSTER = 'RankInCluster'
    LED = 'LED'
 
    ORIGINAL_ID = 'OriginalID'
    
    # Diversity analysis columns
    TOTAL_SEQUENCES = 'TotalSequences'
    TOTAL_READS = 'TotalReads'
    TOTAL_RPU = 'TotalRPU'
    AVERAGE_LED = 'AverageLED'
    SID = 'SID'
    
    # Motif analysis columns
    POPULATION_NUMBER = 'PopulationNumber'
    POPULATION_NAME = 'PopulationName'
    PERCENTAGE = 'Percentage'
    MOTIF = 'Motif'
    ALIAS = 'Alias'
    
    ENRICHMENT = 'Enrichment'
    COMPARISON = 'Comparison'
    QUERY = 'Query'
    
        # motif discovery output
    P_VALUE = "P"
    Z_SCORE = "Z"
    ZZ_SCORE = "ZZ"
    RATIO = "R"
    MOTIF_LENGTH = "Motif_Length"
    SEQ_COUNT = "SeqCount"
    OBSERVED_COUNT = "ObservedCount"
    EXPECTED_COUNT = "ExpectedCount"