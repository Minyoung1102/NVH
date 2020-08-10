clc, clear
%piece = [4 5 6 7 9 11 13]
piece = [9]
for i= piece
    path = ['C:\Users\82108\Desktop\data\Reference_639_grouping_tight\' num2str(i) '\']
    load_img = ['C:\Users\82108\Desktop\data\Reference_639_grouping_tight\' num2str(i) '\*.jpg']

    list=dir(load_img);

    for j = 1:length(list)
        
        file_name = list(j).name
        
        if file_name(1:2) ~= "_p"
            old_img=imread([path file_name]); % 컬러 이미지로 읽어와짐
            
            img= imresize(old_img,[128 128])% binary 이미지로 변환

            %시각화
            f1 = figure
            imshow(img)

            I2 = img(1:64,65:128)

            for m=1:64

                for n= 1:(65-m)*tan((90-360/i)/180*pi)

                  I2(m,n)=255;
                end

            end
            I2 = rot90(I2, 1);
            I2 = fliplr(I2);  
            if i= 9
              
              I2=I2(13:76,1:64);
            end
             if i= 11
              
             I2=I2(36:99,1:64);
            end
             if i= 13
              
             I2=I2(58:121,1:64);
            end
   
              


            %I2 시각화
            f2 = figure
            imshow(I2)

            new_file_name= [path 'p_' file_name]
            imwrite(I2,new_file_name)

            close(f1)
            close(f2)
        end
    end
end