def get_dir(job):
    directory = ''
    count = 1
    for frag in job:
        directory += str(frag)
        if(count < len(job)):
            directory += '_'
        count += 1
    print(directory)
    return(directory)

