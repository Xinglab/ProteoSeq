import os,sys,argparse

def parse_attributes(string):
    d = {}
    for i in string.split(';'):
        if i == '':
            continue
        i_list = i.strip().split(' ')
        if len(i_list) < 2:
            print(string)
            continue
        d.setdefault(i_list[0], []).append(i_list[1].strip('"'))
    return d

def tx_tag_check(tx_tag_list, required_tx_tag_list):
    if not required_tx_tag_list:
        return True
    if len(required_tx_tag_list) == 0:
        return True
    if required_tx_tag_list == ['']:
        return True
    if len(tx_tag_list) == 0:
        return False
    for required_tx_tag in required_tx_tag_list:
        if required_tx_tag not in tx_tag_list: #Some required tag is not found
            return False
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select transcripts from given GTF with certain tags.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file. Required', type=str, required=True)
    parser.add_argument('--output_discarded', help='Output additional GTF file of discarded lines', type=str)
    parser.add_argument('--tx_tag', help='Use CDS of transcripts with selected tags, comma splitted, set off to include all transcripts, default=None', type=str)
    parser.add_argument('--tx_type', help='Use CDS of transcripts with a selected biotype, set off to include all transcripts, default=None', type=str)
    parser.add_argument('--tx_wCDS', help='Kept transcripts with CDS coordinates, default=True', type=str, default='True', choices=['True', 'False'])

    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_gtf = args.output_gtf
    if args.output_discarded:
        output_discarded = args.output_discarded
    else:
        output_discarded = ''

    if (not args.tx_tag) or args.tx_tag == 'False':
        required_tx_tag_list = None
    else:
        required_tx_tag_list = args.tx_tag.split(',')

    if (not args.tx_type) or args.tx_type == 'False':
        required_tx_type = None
    else:
        required_tx_type = args.tx_type

    if args.tx_wCDS == 'True':
        tx_wCDS = True
    else:
        tx_wCDS = False

    output = open(output_gtf, 'w')
    if output_discarded != '':
        output2 = open(output_discarded, 'w')
    else:
        output2 = None

    tx_list_wCDS = set()
    if tx_wCDS:
        for line in open(input_gtf, 'r'):
            if line.startswith('#'):
                continue
            arr = line.rstrip('\n').split('\t')
            if arr[2] == 'CDS':
                d = parse_attributes(arr[8])
                tx_id = d.get('transcript_id', [''])[0]
                if tx_id != '':
                    tx_list_wCDS.add(tx_id)

    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        tx_id = d.get('transcript_id', [''])[0]
        if tx_id == '':
            if output2:
                output2.write(line)
            continue
        if tx_wCDS and (tx_id not in tx_list_wCDS):
            if output2:
                output2.write(line)
            continue
        tx_type = d.get('transcript_type', [''])[0]
        if tx_type == '':
            if output2:
                output2.write(line)
            continue
        if required_tx_type and tx_type != required_tx_type: #Tx does not have required biotype
            if output2:
                output2.write(line)
            continue
        tx_tag_list = d.get('tag', [])
        if not tx_tag_check(tx_tag_list, required_tx_tag_list):
            if output2:
                output2.write(line)
            continue
        output.write(line)

    output.close()
    if output2:
        output2.close()


